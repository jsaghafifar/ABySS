package abyss.lphybeast.tobeast.generator;

import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.inference.operator.kernel.Transform;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import abyss.evolution.substitution.ABySSubstitutionModel;
import abyss.NonReversible;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

public class NonReversibleToBEAST implements GeneratorToBEAST<NonReversible, ABySSubstitutionModel> {
    @Override
    public ABySSubstitutionModel generatorToBEAST(NonReversible nq, BEASTInterface value, BEASTContext context) {

        ABySSubstitutionModel beastNQ = new ABySSubstitutionModel();

        Value<Double[]> rates = nq.getRates();
        Value<Double[]> freqs = nq.getFreq();

        RealParameter ratesParameter = (RealParameter)context.getBEASTObject(rates);
        RealParameter freqParameter = (RealParameter)context.getBEASTObject(freqs);

        int numStates = freqs.value().length;
        char[] states = new char[numStates];

        if (numStates == 4) {
            states = new char[] {'A', 'C', 'G', 'T'};
        } else if (numStates == 20) {
            states = new char[] {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
        } else {
            for (int i=0; i<numStates; i++) {
                states[i] = (char) i;
            }
        }

        String[] keysArray = new String[numStates * numStates - numStates];
        int x = 0;
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                if (j != i) {
                    keysArray[x] = new String(new char[]{states[i], states[j]});
                    x++;
                }
            }
        }
        String keys = String.join(" ", keysArray);
        String stateNames = new String(states).replace("", " ").trim();
        ratesParameter.setInputValue("keys", keys);

        BooleanParameter rateIndicatorParameter = (BooleanParameter)context.getBEASTObject(nq.getIndicators());
        rateIndicatorParameter.setInputValue("keys", keys);
        rateIndicatorParameter.initAndValidate();
        context.addSkipLoggable(rateIndicatorParameter);

        boolean symmetricParameter = false;

        List<Transform> rateTransforms = new ArrayList<>();
        rateTransforms.add(addLogTransform(ratesParameter));
        addAVMNOperator(context, rateTransforms);

        List<Transform> freqTransforms = new ArrayList<>();
        freqTransforms.add(addLogTransform(freqParameter));
        addAVMNOperator(context, freqTransforms);

        ratesParameter.initAndValidate();
        context.addSkipLoggable(ratesParameter);
        context.addSkipLoggable(freqParameter);

        beastNQ.setInputValue("rates", ratesParameter);
        beastNQ.setInputValue("frequencies", BEASTContext.createBEASTFrequencies((RealParameter) context.getBEASTObject(freqs), stateNames));
        beastNQ.setInputValue("rateIndicator", rateIndicatorParameter);
        beastNQ.setInputValue("symmetric", symmetricParameter);
        beastNQ.initAndValidate();
        return beastNQ;
    }

    private void addAVMNOperator(BEASTContext context, List<Transform> transforms) {
        AdaptableVarianceMultivariateNormalOperator operator = new AdaptableVarianceMultivariateNormalOperator();
        operator.initByName("weight", 1.0,"coefficient", 1.0,"scaleFactor", 1.0,"beta", 0.05,
                "initial", 1000,"burnin",500, "allowNonsense", false, "transformations", transforms);
        operator.setID("AVMN"); //TODO fix to differentiate rates/freqs?
        context.addExtraOperator(operator);
    }

    private Transform addLogTransform(RealParameter parameter) {
        List<Function> logFunc = new ArrayList<>();
        logFunc.add(parameter);
        Transform.LogTransform logTransform = new Transform.LogTransform();
        logTransform.initByName("f", logFunc);
        logTransform.setID("logTrans."+parameter.getID());
        return logTransform;
    }

    @Override
    public Class<NonReversible> getGeneratorClass() { return NonReversible.class; }

    @Override
    public Class<ABySSubstitutionModel> getBEASTClass() {
        return ABySSubstitutionModel.class;
    }
}
