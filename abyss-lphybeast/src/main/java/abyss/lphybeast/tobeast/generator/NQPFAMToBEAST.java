package abyss.lphybeast.tobeast.generator;

import beast.base.core.BEASTInterface;
import beast.base.evolution.substitutionmodel.ComplexSubstitutionModel;
import beast.base.inference.parameter.RealParameter;
import abyss.NQPFAM;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

/**
 * @author Jasmine Saghafifar
 */
public class NQPFAMToBEAST implements GeneratorToBEAST<NQPFAM, ComplexSubstitutionModel> {
    @Override
    public ComplexSubstitutionModel generatorToBEAST(NQPFAM nqpfam, BEASTInterface value, BEASTContext context) {

        ComplexSubstitutionModel beastNQPFAM = new ComplexSubstitutionModel();

        RealParameter freqParameter = new RealParameter(nqpfam.getFreq());
        RealParameter ratesParameter = new RealParameter(nqpfam.getRates());

        int numStates = nqpfam.getFreq().length;
        char[] states = new char[] {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
        String[] keysList = new String[numStates * numStates - numStates];
        int x = 0;
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                if (j != i) {
                    keysList[x] = new String(new char[]{states[i], states[j]});
                    x++;
                }
            }
        }
        String keys = String.join(" ", keysList);
        ratesParameter.setInputValue("keys", keys);
        ratesParameter.initAndValidate();

        beastNQPFAM.setInputValue("rates", ratesParameter);
        beastNQPFAM.setInputValue("frequencies", BEASTContext.createBEASTFrequencies(freqParameter, "A C D E F G H I K L M N P Q R S T V W Y"));
        beastNQPFAM.initAndValidate();
        return beastNQPFAM;
    }

    @Override
    public Class<NQPFAM> getGeneratorClass() { return NQPFAM.class; }

    @Override
    public Class<ComplexSubstitutionModel> getBEASTClass() {
        return ComplexSubstitutionModel.class;
    }
}
