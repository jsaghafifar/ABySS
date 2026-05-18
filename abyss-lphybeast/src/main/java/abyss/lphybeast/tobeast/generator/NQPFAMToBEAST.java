package abyss.lphybeast.tobeast.generator;

import abyss.substitutionmodel.ABySSubstitutionModel;
import beast.base.core.BEASTInterface;
import abyss.NQPFAM;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealVectorParam;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import static abyss.lphybeast.tobeast.generator.NonReversibleToBEAST.getRateKeys;
import static abyss.lphybeast.tobeast.generator.NonReversibleToBEAST.getStates;

/**
 * @author Jasmine Saghafifar
 */
public class NQPFAMToBEAST implements GeneratorToBEAST<NQPFAM, ABySSubstitutionModel> {
    @Override
    public ABySSubstitutionModel generatorToBEAST(NQPFAM nqpfam, BEASTInterface value, BEASTContext context) {

        ABySSubstitutionModel beastNQPFAM = new ABySSubstitutionModel();

        RealVectorParam<PositiveReal> ratesParameter = new RealVectorParam<>(nqpfam.getRates(), PositiveReal.INSTANCE);
        int numStates = nqpfam.getNumStates();

        char[] states = getStates(context, numStates);
        String keys = getRateKeys(states, numStates, false);
        ratesParameter.setInputValue("keys", keys);
        ratesParameter.initAndValidate();

        beastNQPFAM.setInputValue("rates", ratesParameter);
        beastNQPFAM.initAndValidate();
        return beastNQPFAM;
    }

    @Override
    public Class<NQPFAM> getGeneratorClass() { return NQPFAM.class; }

    @Override
    public Class<ABySSubstitutionModel> getBEASTClass() {
        return ABySSubstitutionModel.class;
    }
}
