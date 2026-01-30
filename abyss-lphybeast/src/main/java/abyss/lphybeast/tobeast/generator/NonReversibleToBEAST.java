package abyss.lphybeast.tobeast.generator;

import abyss.MixedAlignment;
import abyss.distributions.EigenFriendlyQPrior;
import abyss.distributions.SVSPrior;
import abyss.logger.ABySSFrequencyLogger;
import abyss.logger.DetailedBalanceLogger;
import abyss.logger.NetFluxLogger;
import abyss.logger.RootMeanSquareLogger;
import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.inference.operator.BitFlipOperator;
import beast.base.inference.operator.kernel.Transform;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import abyss.substitutionmodel.ABySSubstitutionModel;
import abyss.NonReversible;
import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.likelihood.PhyloCTMC;
import lphy.core.model.Value;
import lphy.core.vectorization.array.BooleanArray;
import lphy.core.vectorization.array.DoubleArray;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Jasmine Saghafifar
 */

public class NonReversibleToBEAST implements GeneratorToBEAST<NonReversible, ABySSubstitutionModel> {
    @Override
    public ABySSubstitutionModel generatorToBEAST(NonReversible nq, BEASTInterface value, BEASTContext context) {
        ABySSubstitutionModel beastNQ = new ABySSubstitutionModel();
        int numStates;
        Value<Boolean> symmetric = nq.getSymmetric();
        Value<Double[]> rates = nq.getRates();
        RealParameter ratesParameter = (RealParameter)context.getBEASTObject(rates);

        // get number of states for keys
        if (symmetric.value()) {
            numStates = nq.getFreq().value().length;
        } else {
            double root = (-1 - Math.sqrt(1+4* rates.value().length))/2;
            numStates = (int) Math.abs(root);
        }

        // make keys for freq and rates
        char[] states = getStates(context, numStates);
        String stateNames = new String(states).replace("", " ").trim();
        String keys = getRateKeys(states, numStates, symmetric.value());

        // init reversible Q freqs
        if (symmetric.value()) {
            RealParameter freqParameter = (RealParameter)context.getBEASTObject(nq.getFreq());
            beastNQ.setInputValue("frequencies", BEASTContext.createBEASTFrequencies(freqParameter, stateNames));
        }

        ratesParameter.setInputValue("keys", keys);
        ratesParameter.initAndValidate();
        ratesParameter.setID(rates.getId());

        // add rate indicators and non-uniform operator
        if (nq.getIndicators() != null) {
            Value<Boolean[]> indicators = nq.getIndicators();

            BEASTInterface bi = context.getBEASTObject(indicators);
            BooleanParameter rateIndicatorParameter;
            if (bi instanceof SVSPrior prior) {
                rateIndicatorParameter = prior.indicatorsInput.get();
            } else if (bi instanceof BooleanParameter) {
                rateIndicatorParameter = (BooleanParameter) bi;
            } else throw new UnsupportedOperationException();

            rateIndicatorParameter.setInputValue("keys", keys);
            rateIndicatorParameter.initAndValidate();
            rateIndicatorParameter.setID(indicators.getId());
            if (!symmetric.value()) createEigenFriendlyQPrior(context, ratesParameter, rateIndicatorParameter, numStates, value.getID());

            if (!(indicators.getGenerator() instanceof BooleanArray | indicators.getGenerator() instanceof ConnectedSVS))
                addNonUniformBitFlipOperator(context, rateIndicatorParameter, 2.0);

            beastNQ.setInputValue("rateIndicator", rateIndicatorParameter);
        } else if (!symmetric.value())
            createEigenFriendlyQPrior(context, ratesParameter, null, numStates, value.getID());

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

    // public static methods
    public static char[] getStates(BEASTContext context, int numStates) {
        char[] states = new char[numStates];

        if (((context.getAlignments().get(0).getGenerator() instanceof PhyloCTMC phyloCTMC &&
                phyloCTMC.getDataType() == SequenceType.NUCLEOTIDE) ||
                context.getAlignments().get(0).getGenerator() instanceof MixedAlignment mixedAlignment &&
                        mixedAlignment.getDataType() == SequenceType.NUCLEOTIDE)
                && numStates == 4)
            states = new char[] {'A', 'C', 'G', 'T'};
        else if (((context.getAlignments().get(0).getGenerator() instanceof PhyloCTMC phyloCTMC &&
                phyloCTMC.getDataType() == SequenceType.AMINO_ACID) ||
                context.getAlignments().get(0).getGenerator() instanceof MixedAlignment mixedAlignment &&
                        mixedAlignment.getDataType() == SequenceType.AMINO_ACID)
                && numStates == 20)
            states = new char[] {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
        else {
            for (int i=0; i<numStates; i++) {
                states[i] = Character.forDigit(i, 10);
            }
        }
        return states;
    }

    public static String getRateKeys(char[] states, int numStates, Boolean symmetric) {
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

    public static void addNonUniformBitFlipOperator(BEASTContext context, BooleanParameter parameter, double weight) {
        context.addSkipOperator(parameter);
        BitFlipOperator operator = new BitFlipOperator();
        operator.setID(parameter.getID() + ".bitFlip");
        operator.initByName("weight", weight, "parameter", parameter, "uniform", true);
        operator.setID(parameter.getID() + ".bitFlip");
        context.addExtraOperator(operator);
    }

    // private operator methods
    private void addAVMNOperator(BEASTContext context, List<Transform> transforms, double weight, String id) {
        AdaptableVarianceMultivariateNormalOperator operator = new AdaptableVarianceMultivariateNormalOperator();
        operator.initByName("weight", weight,"coefficient", 1.0,"scaleFactor", 1.0,"beta", 0.05,
                "initial", 1000,"burnin",500, "allowNonsense", false, "transformations", transforms);
        operator.setID("AVMN."+id);
        context.addExtraOperator(operator);
    }

    private Transform addLogConstrainedSumTransform(RealParameter parameter) {
        List<Function> logFunc = new ArrayList<>();
        logFunc.add(parameter);
        Transform.LogConstrainedSumTransform logConstrainedSumTransform = new Transform.LogConstrainedSumTransform();
        logConstrainedSumTransform.initByName("f", logFunc);
        logConstrainedSumTransform.setID("logSumTrans."+parameter.getID());
        return logConstrainedSumTransform;
    }

    // private logger methods
    private void addFreqLogger(BEASTContext context, ABySSubstitutionModel model, String keys, String id) {
        ABySSFrequencyLogger frequencyLogger = new ABySSFrequencyLogger();
        frequencyLogger.setInputValue("model", model);
        frequencyLogger.setInputValue("keys", keys);
        frequencyLogger.initAndValidate();
        if (context.getAlignments().get(0).getGenerator() instanceof MixedAlignment)
        frequencyLogger.setID(id+".equilibrium.freq"); else frequencyLogger.setID("equilibrium.freq");

        context.addExtraLoggable(frequencyLogger);
    }

    private void addNonRevLoggers(BEASTContext context, ABySSubstitutionModel model, String keys, String states, Integer numStates) {
        DetailedBalanceLogger balancesLogger = new DetailedBalanceLogger();
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
            NetFluxLogger fluxesLogger = new NetFluxLogger();
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

    // private prior methods
    private void createEigenFriendlyQPrior(BEASTContext context, Function rates, BooleanParameter indicators, Integer numStates, String id) {
        // ensure that nonreversible Q matrix is eigen-friendly
        EigenFriendlyQPrior qPrior = new EigenFriendlyQPrior();
        qPrior.setInputValue("rates", rates);
        if (indicators != null) {
            qPrior.setInputValue("rateIndicator", indicators);
        }
        qPrior.setInputValue("nrOfStates", numStates);
        qPrior.initAndValidate();
        qPrior.setID("eigenfriendly"+id+".prior");
        context.addBEASTObject(qPrior, null);

    }

    @Override
    public Class<NonReversible> getGeneratorClass() { return NonReversible.class; }

    @Override
    public Class<ABySSubstitutionModel> getBEASTClass() {
        return ABySSubstitutionModel.class;
    }
}
