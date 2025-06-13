package abyss.distributions;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.substitutionmodel.ColtEigenSystem;
import beast.base.evolution.substitutionmodel.EigenSystem;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import jdk.jfr.Description;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static abyss.evolution.substitution.ABySSubstitutionModel.setupUnnormNonrevQ;

@Description("Check can converge, reject proposed Q matrix otherwise")
public class EigenFriendlyQPrior extends Distribution {

    final public Input<Function> ratesInput = new Input<>("rates",
            "rates parameter", Input.Validate.REQUIRED);
    final public Input<BooleanParameter> indicatorsInput = new Input<>("rateIndicator",
            "rates to indicate the presence or absence of transition matrix entries", Input.Validate.OPTIONAL);
    final public Input<Integer> nrOfStatesInput = new Input<>("nrOfStates",
            "number of states parameter", Input.Validate.REQUIRED);
    protected Function rates;
    protected double[][] Qm;
    protected Integer nrOfStates;
    protected EigenSystem eigenSystem;


    @Override
    public double calculateLogP() {
        logP = 0.0;
        Function rates = this.ratesInput.get();
        nrOfStates = this.nrOfStatesInput.get();
        double[] relativeRates = new double[rates.getDimension()];
        if (this.indicatorsInput.get() != null) {
            for (int i = 0; i < relativeRates.length; i++) {
                relativeRates[i] = rates.getArrayValue(i) * (indicatorsInput.get().getValue(i)?1.:0.);
            }
        } else for (int i = 0; i < relativeRates.length; i++) {
            relativeRates[i] = rates.getArrayValue(i);
        }

        Qm = setupUnnormNonrevQ(relativeRates, nrOfStates);
        eigenSystem = new ColtEigenSystem(nrOfStates);
        try { eigenSystem.decomposeMatrix(Qm); } catch (Exception exception) {
            logP = Double.NEGATIVE_INFINITY;
        }
        return logP;
    }

    @Override
    public List<String> getArguments() {
        List<String> args = new ArrayList<>();
        args.add(((RealParameter)ratesInput.get()).getID());
        if (indicatorsInput.get() != null) {
            args.add((indicatorsInput.get().getID()));
        }
        return args;
    }

    @Override
    public List<String> getConditions() {
        List<String> conds = new ArrayList<>();
        return conds;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void initAndValidate() {
        if (indicatorsInput.get() != null &&
                indicatorsInput.get().getDimension() != ratesInput.get().getDimension()) {
            throw new RuntimeException("Indicators must be same size as rates parameter but it was dimension " +
                    indicatorsInput.get().getDimension());
        }
    }
}
