package abyss.substitutionmodel;

import beast.base.core.*;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;

/**
 * @author Jasmine Saghafifar
 */

@Description("Equilibrium frequency logger")
public class ABySSFrequencyLogger extends CalculationNode implements Loggable, Function {
    final public Input<ABySSubstitutionModel> modelInput;
    final public Input<String[]> keysInput;
    protected ABySSubstitutionModel model;
    protected String[] keys;
    private static final double DEFAULT_BRANCH_LENGTH = 100000;

    public ABySSFrequencyLogger() {
        this.modelInput = new Input("model", "ABYSS SVS general substitution model.", Input.Validate.REQUIRED);
        this.keysInput = new Input("keys", "freq state names");
    }

    public void initAndValidate() {
        this.model = this.modelInput.get();
        this.keys = this.keysInput.get();
    }

    public void init(PrintStream out) {
        // substModel + freq + state name of j = param name
        String id = "freq.";
        String stateName;

        for (int i = 0; i < this.model.getStateCount(); i++) {
            stateName = this.keys != null ? keys[i] : String.valueOf(i);
            out.print(id + stateName + "\t");
        }

    }

    public void log(long sample, PrintStream out) {
        double[] p = getDoubleValues();
        for (int i = 0; i < this.model.getStateCount(); i++) {
            out.print(p[i] + "\t");
        }

    }

    private double[] getFrequencies(ABySSubstitutionModel model) {
        int nrOfStates = model.getStateCount();
        double t = DEFAULT_BRANCH_LENGTH;
        boolean equilibrium = false;
        double[] p;
        double[] f = new double[nrOfStates];

        while (!equilibrium) {
            p = new double[nrOfStates * nrOfStates];
            model.getTransitionProbabilities(null, 1.0, 0., t, p);
            boolean reached = true;
            for (int i = 0; i < nrOfStates; i++) {
                f[i] = p[i];
                for (int j = 1; j < nrOfStates; j++) {
                    if (p[i] - p[j * nrOfStates + i] > 1e-6) {
                        reached = false;
                        break;
                    }
                }
                if (!reached) break;
            }
            if (reached) equilibrium = true;
            t *= 10;
        }
        return f;
    }

    public void close(PrintStream close) {
    }

    @Override
    public int getDimension() {
        return this.model.getStateCount();
    }

    @Override
    public double[] getDoubleValues() {
        return getFrequencies(this.model);
    }

    @Override
    public double getArrayValue(int dim) {
        return getDoubleValues()[dim];
    }
}
