package abyss.distributions;

import beast.base.core.Input;
import beast.base.evolution.substitutionmodel.ColtEigenSystem;
import beast.base.evolution.substitutionmodel.EigenSystem;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.BoolVectorParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealVector;
import jdk.jfr.Description;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static abyss.substitutionmodel.ABySSubstitutionModel.setupUnnormNonrevQ;

/**
 * @author Jasmine Saghafifar
 */

@Description("Check can converge, reject proposed Q matrix otherwise")
public class EigenFriendlyQPrior extends Distribution {

    final public Input<RealVectorParam<NonNegativeReal>> ratesInput = new Input<>("rates", // nonnegative?
            "rates parameter", Input.Validate.REQUIRED);
    final public Input<BoolVectorParam> indicatorsInput = new Input<>("rateIndicator",
            "rates to indicate the presence or absence of transition matrix entries", Input.Validate.OPTIONAL);
    final public Input<Integer> nrOfStatesInput = new Input<>("nrOfStates",
            "number of states parameter", Input.Validate.REQUIRED);
    protected RealVector<NonNegativeReal> rates;
    protected double[][] Qm;
    protected Integer nrOfStates;
    protected EigenSystem eigenSystem;


    @Override
    public double calculateLogP() {
        logP = 0.0;
        RealVector<NonNegativeReal> rates = this.ratesInput.get();
        nrOfStates = this.nrOfStatesInput.get();
        double[] relativeRates = new double[rates.size()];
        if (this.indicatorsInput.get() != null) {
            for (int i = 0; i < relativeRates.length; i++) {
                relativeRates[i] = rates.get(i) * (indicatorsInput.get().get(i)?1.:0.);
            }
        } else for (int i = 0; i < relativeRates.length; i++) {
            relativeRates[i] = rates.get(i);
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
        args.add((ratesInput.get()).getID());
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
                indicatorsInput.get().size() != ratesInput.get().size()) {
            throw new RuntimeException("Indicators must be same size as rates parameter but it was dimension " +
                    indicatorsInput.get().size());
        }
    }
}
