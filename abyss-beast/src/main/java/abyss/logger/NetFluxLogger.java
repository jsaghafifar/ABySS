package abyss.logger;

import abyss.substitutionmodel.ABySSubstitutionModel;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Jasmine Saghafifar
 */

@Description("Logger for net fluxes in nonreversible Q")
public class NetFluxLogger extends CalculationNode implements Loggable, Function {
    final public Input<ABySSubstitutionModel> modelInput;
    final public Input<String> stateNamesInput;
    protected ABySSubstitutionModel model;
    protected List<String> stateNames;

    public NetFluxLogger() {
        this.modelInput = new Input<>("model", "ABYSS SVS general substitution model.",
                Input.Validate.REQUIRED);
        this.stateNamesInput = new Input<>("states", "State names");
    }

    public void initAndValidate() {
        this.model = this.modelInput.get();
        if (model.getStateCount() != 4)
            throw new IllegalArgumentException("Logger currently only supports substitution with 4 states.");

        if (this.stateNamesInput.get() != null) {
            String[] keysArr = stateNamesInput.get().split(" ");
            List<String> stateNames = Collections.unmodifiableList(Arrays.asList(keysArr));
            if (this.stateNames.size() != getDimension())
                throw new IllegalArgumentException("Number of keys must match nrOfStates for model. " +
                        "Model nrOfStates = " + model.getStateCount() + ", but stateNames.size() = " +
                        this.stateNames.size());
            this.stateNames = stateNames;
        }
    }

    public void init(PrintStream out) {
        // logger id + cycle name = param name
        String id = getID();
        
        String[] states = this.stateNames != null ? stateNames.toArray(new String[0]) : new String[]{"1","2","3","4"};
        int[][] cycles = getCycles();

        for (int[] cycle : cycles) {
            String cycleName = "";
            for (int j = 0; j < cycle.length; j++) {
                String step = states[cycle[j]];
                cycleName = cycleName.concat(step);
            }

            out.print(id + '.' + cycleName + "\t");
        }

    }

    public void log(long sample, PrintStream out) {
        double[] p = getDoubleValues();
        for (int i = 0; i < this.model.getStateCount(); i++) {
            out.print(p[i] + "\t");
        }

    }

    private double[] getNetFluxes() {
        model.setupRateMatrix();
        double[][] Q = model.getRateMatrix();
        int[][] cycles = getCycles();

        double[] f = new double[cycles.length];
        int x = 0;

        // based on the detailed balance property (Q_ij*pi_i=Q_ji*pi_j),
        // cycles should have zero flux when time-reversible
        // e.g. A-C-G-T, forward Q_AC*Q_CG*Q_GT*Q_TA = backward Q_AT*Q_TG*Q_GC*Q_CA
        // presence of net flux indicates time-nonreversibility
        for (int[] cycle : cycles) {
            double forward = 1;
            double backward = 1;
            for (int step = 0; step < cycle.length; step++) {
                int i = cycle[step];
                int j = step < cycle.length - 1 ? cycle[step+1] : cycle[1];
                forward *= Q[i][j];
                backward *= Q[j][i];
            }
            f[x] = forward - backward;
        }

        return f;
    }

    private int[][] getCycles() { // TODO make generalisable for different nrOfStates (low priority)
        int[][] cycles = new int[][]{ {1,2,3},  {1,2,4},  {1,3,4},  {2,3,4},
                                      {1,2,3,4},{1,2,4,3},{1,3,2,4}        };
        return cycles;
    }

    public void close(PrintStream close) {
    }

    @Override
    public int getDimension() {
        return 7;
    }

    @Override
    public double[] getDoubleValues() {
        return getNetFluxes();
    }

    @Override
    public double getArrayValue(int dim) {
        return getDoubleValues()[dim];
    }
}
