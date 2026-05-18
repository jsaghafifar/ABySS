package abyss.lphybeast.tobeast.generator;

import abyss.ClassicNonReversible;
import abyss.distributions.SVSPrior;
import abyss.logger.*;
import beast.base.core.BEASTInterface;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.inference.operator.kernel.Transform;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.BoolVectorParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel;
import lphy.core.model.Value;
import lphy.core.vectorization.array.DoubleArray;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static abyss.lphybeast.tobeast.generator.NonReversibleToBEAST.getStates;
import static abyss.lphybeast.tobeast.generator.NonReversibleToBEAST.getRateKeys;

/**
 * @author Jasmine Saghafifar
 */

public class ClassicNonReversibleToBEAST implements GeneratorToBEAST<ClassicNonReversible, SVSGeneralSubstitutionModel> {
    @Override
    public SVSGeneralSubstitutionModel generatorToBEAST(ClassicNonReversible nq, BEASTInterface value, BEASTContext context) {
        SVSGeneralSubstitutionModel beastNQ = new SVSGeneralSubstitutionModel();
        int numStates;
        Value<Boolean> symmetric = nq.getSymmetric();
        Value<Double[]> rates = nq.getRates();
        RealVectorParam<PositiveReal> ratesParameter = (RealVectorParam<PositiveReal>) context.getAsRealVector(rates);

        // get number of states for keys
        numStates = nq.getFreq().value().length;

        // make keys for freq and rates
        char[] states = getStates(context, numStates);
        String stateNames = new String(states).replace("", " ").trim();
        String keys = getRateKeys(states, numStates, symmetric.value());

        // init Q freqs
        SimplexParam freqParameter = (SimplexParam) context.getBEASTObject(nq.getFreq());
        freqParameter.setInputValue("keys", stateNames);
        freqParameter.initAndValidate();
        freqParameter.setID(nq.getFreq().getId());
        beastNQ.setInputValue("frequencies", freqParameter);

        ratesParameter.setInputValue("keys", keys);
        ratesParameter.initAndValidate();
        ratesParameter.setID(rates.getId());

        // add rate indicators and non-uniform operator
        if (nq.getIndicators() != null) {
            Value<Boolean[]> indicators = nq.getIndicators();

            BEASTInterface bi = context.getBEASTObject(indicators);
            BoolVectorParam rateIndicatorParameter;
            if (bi instanceof SVSPrior prior) {
                rateIndicatorParameter = prior.indicatorsInput.get();
            } else if (bi instanceof BoolVectorParam) {
                rateIndicatorParameter = (BoolVectorParam) bi;
            } else throw new UnsupportedOperationException();

            rateIndicatorParameter = (BoolVectorParam) context.getBEASTObject(indicators);
            rateIndicatorParameter.setInputValue("keys", keys);
            rateIndicatorParameter.initAndValidate();
            rateIndicatorParameter.setID(indicators.getId());
            beastNQ.setInputValue("rateIndicator", rateIndicatorParameter);
        }

        if (!(nq.getRates().getGenerator() instanceof DoubleArray)) {
            List<Transform> rateTransforms = new ArrayList<>();
            rateTransforms.add(addLogConstrainedSumTransform(ratesParameter));
            addAVMNOperator(context, rateTransforms, 10.0, ratesParameter.getID());
        }

        beastNQ.setInputValue("rates", ratesParameter);
        beastNQ.setInputValue("symmetric", symmetric.value());
        beastNQ.initAndValidate();
        beastNQ.setID(value.getID());
        if (!symmetric.value()) addFreqLogger(context, beastNQ, stateNames, value.getID());
        addNonRevLoggers(context, beastNQ, getRateKeys(states, numStates, true), stateNames, numStates);

        return beastNQ;
    }

    private void addAVMNOperator(BEASTContext context, List<Transform> transforms, double weight, String id) {
        AdaptableVarianceMultivariateNormalOperator operator = new AdaptableVarianceMultivariateNormalOperator();
        operator.initByName("weight", weight,"coefficient", 1.0,"scaleFactor", 1.0,"beta", 0.05,
                "initial", 1000,"burnin",500, "allowNonsense", false, "transformations", transforms);
        operator.setID("AVMN."+id);
        context.addExtraOperator(operator);
    }

    private Transform addLogConstrainedSumTransform(RealVectorParam<?> parameter) {
        List<RealVectorParam<?>> logFunc = new ArrayList<>();
        logFunc.add(parameter);
        Transform.LogConstrainedSumTransform logConstrainedSumTransform = new Transform.LogConstrainedSumTransform();
        logConstrainedSumTransform.initByName("f", logFunc);
        logConstrainedSumTransform.setID("logSumTrans."+parameter.getID());
        return logConstrainedSumTransform;
    }

    private void addFreqLogger(BEASTContext context, SVSGeneralSubstitutionModel model, String keys, String id) {
        ClassicFrequencyLogger frequencyLogger = new ClassicFrequencyLogger();
        frequencyLogger.setInputValue("model", model);
        frequencyLogger.setInputValue("keys", keys);
        frequencyLogger.initAndValidate();
        frequencyLogger.setID("equilibrium.freq");

        context.addExtraLoggable(frequencyLogger);
    }

    private void addNonRevLoggers(BEASTContext context, SVSGeneralSubstitutionModel model, String keys, String states, Integer numStates) {
        ClassicDetailedBalanceLogger balancesLogger = new ClassicDetailedBalanceLogger();
        balancesLogger.setInputValue("model", model);
        balancesLogger.setInputValue("keys", keys);
        balancesLogger.initAndValidate();
        balancesLogger.setID("balances");
        RootMeanSquareLogger balanceLogger = new RootMeanSquareLogger();
        balanceLogger.setInputValue("logger", balancesLogger);
        balanceLogger.initAndValidate();
        balanceLogger.setID("balance");

        context.addExtraLoggable(balancesLogger);
        context.addExtraLoggable(balanceLogger);

        if (numStates==4) {
            ClassicNetFluxLogger fluxesLogger = new ClassicNetFluxLogger();
            fluxesLogger.setInputValue("model", model);
            fluxesLogger.setInputValue("states", states);
            fluxesLogger.initAndValidate();
            fluxesLogger.setID("fluxes");
            RootMeanSquareLogger fluxLogger = new RootMeanSquareLogger();
            fluxLogger.setInputValue("logger", balancesLogger);
            fluxLogger.initAndValidate();
            fluxLogger.setID("flux");

            context.addExtraLoggable(fluxesLogger);
            context.addExtraLoggable(fluxLogger);
        }

    }

    @Override
    public Class<ClassicNonReversible> getGeneratorClass() { return ClassicNonReversible.class; }

    @Override
    public Class<SVSGeneralSubstitutionModel> getBEASTClass() { return SVSGeneralSubstitutionModel.class; }
}
