package abyss.lphybeast.tobeast.generator;

import abyss.evolution.substitution.ABySSFrequencyLogger;
import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.inference.operator.kernel.Transform;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import abyss.evolution.substitution.ABySSubstitutionModel;
import abyss.NonReversible;
import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.likelihood.PhyloCTMC;
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
        Value<Boolean[]> indicators = nq.getIndicators();
        Value<Boolean> symmetric = nq.getSymmetric();

        int numStates = freqs.value().length;

        RealParameter ratesParameter = (RealParameter)context.getBEASTObject(rates);
        RealParameter freqParameter = (RealParameter)context.getBEASTObject(freqs);
        BooleanParameter rateIndicatorParameter = (BooleanParameter)context.getBEASTObject(indicators);

        char[] states = getStateNames(context, numStates);
        String stateNames = new String(states).replace("", " ").trim();
        String keys = getRateKeys(states, numStates, symmetric.value());

        String[] statesList = new String[states.length];
        for (int i = 0; i < states.length; i++) {
            statesList[i] = String.valueOf(states[i]);
        }

        ratesParameter.setInputValue("keys", keys);
        ratesParameter.initAndValidate();
        rateIndicatorParameter.setInputValue("keys", keys);
        rateIndicatorParameter.initAndValidate();
        //TODO add bitflip operator spec: uniform=false

        List<Transform> rateTransforms = new ArrayList<>();
        rateTransforms.add(addLogConstrainedSumTransform(ratesParameter));
        addAVMNOperator(context, rateTransforms, 10.0, "rates");

        beastNQ.setInputValue("rates", ratesParameter);
        beastNQ.setInputValue("frequencies", BEASTContext.createBEASTFrequencies(freqParameter, stateNames));
        beastNQ.setInputValue("rateIndicator", rateIndicatorParameter);
        beastNQ.setInputValue("symmetric", symmetric.value());
        beastNQ.initAndValidate();
        addFreqLogger(context, beastNQ, statesList);

        return beastNQ;
    }

    private void addAVMNOperator(BEASTContext context, List<Transform> transforms, double weight, String id) {
        AdaptableVarianceMultivariateNormalOperator operator = new AdaptableVarianceMultivariateNormalOperator();
        operator.initByName("weight", weight,"coefficient", 1.0,"scaleFactor", 1.0,"beta", 0.05,
                "initial", 1000,"burnin",500, "allowNonsense", false, "transformations", transforms);
        operator.setID("AVMN."+id);
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

    private Transform addLogConstrainedSumTransform(RealParameter parameter) {
        List<Function> logFunc = new ArrayList<>();
        logFunc.add(parameter);
        Transform.LogConstrainedSumTransform logConstrainedSumTransform = new Transform.LogConstrainedSumTransform();
        logConstrainedSumTransform.initByName("f", logFunc);
        logConstrainedSumTransform.setID("logSumTrans."+parameter.getID());
        return logConstrainedSumTransform;
    }

    private void addFreqLogger(BEASTContext context, ABySSubstitutionModel model, String[] keys) {
        ABySSFrequencyLogger frequencyLogger = new ABySSFrequencyLogger();
        frequencyLogger.setInputValue("model", model);
//        frequencyLogger.setInputValue("keys", keys);
        frequencyLogger.initAndValidate();

        context.addExtraLoggable(frequencyLogger);
    }

    private char[] getStateNames(BEASTContext context, int numStates) {
        char[] states = new char[numStates];

        if (numStates == 4 && context.getAlignments().get(0).getGenerator() instanceof
                PhyloCTMC phyloCTMC && phyloCTMC.getDataType() == SequenceType.NUCLEOTIDE) {
            states = new char[] {'A', 'C', 'G', 'T'};
        } else if (numStates == 20 && context.getAlignments().get(0).getGenerator() instanceof
                PhyloCTMC phyloCTMC && phyloCTMC.getDataType() == SequenceType.AMINO_ACID) {
            states = new char[] {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
        } else {
            for (int i=0; i<numStates; i++) {
                states[i] = (char) i;
            }
        }
        return states;
    }

    private String getRateKeys(char[] states, int numStates, Boolean symmetric) {
        String[] keysArray;
        int x = 0;
        if (!symmetric) {
            keysArray = new String[numStates * numStates - numStates];
            for (int i = 0; i < numStates; i++) {
                for (int j = 0; j < numStates; j++) {
                    if (j != i) {
                        keysArray[x] = new String(new char[]{states[i], '.', states[j]});
                        x++;
                    }
                }
            }
        } else {
            keysArray = new String[(numStates * numStates - numStates)/2];
            for (int i = 0; i < numStates; i++) {
                for (int j = i + 1; j < numStates; j++) {
                    keysArray[x] = new String(new char[]{states[i], '.', states[j]});
                    x++;
                }
            }
        }

        return String.join(" ", keysArray);
    }

    @Override
    public Class<NonReversible> getGeneratorClass() { return NonReversible.class; }

    @Override
    public Class<ABySSubstitutionModel> getBEASTClass() {
        return ABySSubstitutionModel.class;
    }
}
