package abyss.evolution.substitution;

import beast.base.core.*;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;

@Description("Equilibrium frequency logger")
public class ABySSFrequencyLogger extends CalculationNode implements Loggable, Function {
    final public Input<ABySSubstitutionModel> modelInput;
    final public Input<String[]> keysInput;
    protected ABySSubstitutionModel model;
    protected String[] keys;

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
        int numStates = model.getStateCount();
        double[] p = new double[numStates*numStates];
        model.getTransitionProbabilities(null, 1.0, 0., 1000000, p);
        for (int i = 0; i+1 < numStates; i++) {
            for (int j = 0; j < numStates; j++) {
                if (p[i*numStates + j] - p[(i+1)*numStates + j] > 1e-12) throw new ArithmeticException("Did not reach equilibrium");
            }
        }
        return p;
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
