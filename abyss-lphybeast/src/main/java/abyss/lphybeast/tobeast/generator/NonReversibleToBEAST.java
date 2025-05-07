package abyss.lphybeast.tobeast.generator;

import beast.base.core.BEASTInterface;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import abyss.evolution.substitution.ABySSubstitutionModel;
import abyss.NonReversible;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

public class NonReversibleToBEAST implements GeneratorToBEAST<NonReversible, ABySSubstitutionModel> {
    @Override
    public ABySSubstitutionModel generatorToBEAST(NonReversible nq, BEASTInterface value, BEASTContext context) {

        ABySSubstitutionModel beastNQ = new ABySSubstitutionModel();

//        RealParameter freqParameter = new RealParameter(nq.getFreq().value());
        Value<Double[]> rates = nq.getRates();

        RealParameter ratesParameter = (RealParameter)context.getBEASTObject(rates);

        int numStates = nq.getFreq().value().length;
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
        ratesParameter.initAndValidate();

        BooleanParameter rateIndicatorParameter = (BooleanParameter)context.getBEASTObject(nq.getIndicators());
        rateIndicatorParameter.setInputValue("keys", keys);
        rateIndicatorParameter.initAndValidate();

        boolean symmetricParameter = false;

        beastNQ.setInputValue("rates", ratesParameter);
        beastNQ.setInputValue("frequencies", BEASTContext.createBEASTFrequencies((RealParameter) context.getBEASTObject(nq.getFreq()), stateNames));
        beastNQ.setInputValue("rateIndicator", rateIndicatorParameter);
        beastNQ.setInputValue("symmetric", symmetricParameter);
        beastNQ.initAndValidate();
        return beastNQ;
    }

    @Override
    public Class<NonReversible> getGeneratorClass() { return NonReversible.class; }

    @Override
    public Class<ABySSubstitutionModel> getBEASTClass() {
        return ABySSubstitutionModel.class;
    }
}
