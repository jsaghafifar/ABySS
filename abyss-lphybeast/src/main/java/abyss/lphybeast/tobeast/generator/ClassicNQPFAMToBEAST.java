package abyss.lphybeast.tobeast.generator;

import abyss.ClassicNQPFAM;
import beast.base.core.BEASTInterface;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.BoolVectorParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.Arrays;

import static abyss.lphybeast.tobeast.generator.NonReversibleToBEAST.getRateKeys;
import static abyss.lphybeast.tobeast.generator.NonReversibleToBEAST.getStates;

/**
 * @author Jasmine Saghafifar
 */
public class ClassicNQPFAMToBEAST implements GeneratorToBEAST<ClassicNQPFAM, SVSGeneralSubstitutionModel> {
    @Override
    public SVSGeneralSubstitutionModel generatorToBEAST(ClassicNQPFAM nqpfam, BEASTInterface value, BEASTContext context) {

        SVSGeneralSubstitutionModel beastNQPFAM = new SVSGeneralSubstitutionModel();

        RealVectorParam<PositiveReal> ratesParameter = new RealVectorParam<>(nqpfam.getRates(), PositiveReal.INSTANCE);
        SimplexParam freqParameter = new SimplexParam(nqpfam.getFreq());

        boolean[] b = new boolean[380];
        Arrays.fill(b, true);
        BoolVectorParam rateIndicatorParameter = new BoolVectorParam(b);

        int numStates = nqpfam.getNumStates();
        char[] states = getStates(context, numStates);
        String stateNames = new String(states).replace("", " ").trim();
        String keys = getRateKeys(states, numStates, false);
        ratesParameter.setInputValue("keys", keys);
        ratesParameter.initAndValidate();
        freqParameter.setInputValue("keys", stateNames);
        freqParameter.initAndValidate();

        beastNQPFAM.setInputValue("rates", ratesParameter);
        beastNQPFAM.setInputValue("frequencies", freqParameter);
        beastNQPFAM.setInputValue("rateIndicator",rateIndicatorParameter);
        beastNQPFAM.setInputValue("symmetric",false);

        beastNQPFAM.initAndValidate();
        return beastNQPFAM;
    }

    @Override
    public Class<ClassicNQPFAM> getGeneratorClass() { return ClassicNQPFAM.class; }

    @Override
    public Class<SVSGeneralSubstitutionModel> getBEASTClass() {
        return SVSGeneralSubstitutionModel.class;
    }
}
