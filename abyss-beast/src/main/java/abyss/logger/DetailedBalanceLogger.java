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

@Description("Logger for deviations from detailed balance in nonreversible Q")
public class DetailedBalanceLogger extends CalculationNode implements Loggable, Function {
    final public Input<ABySSubstitutionModel> modelInput;
    final public Input<String> keysInput;
    protected ABySSubstitutionModel model;
    protected List<String> keys;
    protected double[] equilibriumFreqs;

    public DetailedBalanceLogger() {
        this.modelInput = new Input<>("model", "ABYSS SVS general substitution model.",
                Input.Validate.REQUIRED);
        this.keysInput = new Input<>("keys", "Rate names");
    }

    public void initAndValidate() {
        this.model = this.modelInput.get();
        this.equilibriumFreqs = model.getFrequencies();

        if (this.keysInput.get() != null) {
            String[] keysArr = keysInput.get().split(" ");
            List<String> keys = Collections.unmodifiableList(Arrays.asList(keysArr));
            if (keys.size() != getDimension())
                throw new IllegalArgumentException("Number of keys must match dimension. " +
                        "Dimension = " + getDimension() + ", but keys.size() = " + keys.size());
            this.keys = keys;
        }
    }

    public void init(PrintStream out) {
        // logger id + rate name = param name
        String id = getID();
        String rateName;
        for (int i = 0; i < getDimension(); i++) {
            rateName = this.keys != null ? keys.get(i) : String.valueOf(i);
            out.print(id + '.' + rateName + "\t");
        }

    }

    public void log(long sample, PrintStream out) {
        double[] p = getDoubleValues();
        for (int i = 0; i < getDimension(); i++) {
            out.print(p[i] + "\t");
        }

    }

    private double[] getDeviations() {
        int nrOfStates = model.getStateCount();
        double[] d = new double[(nrOfStates*nrOfStates-nrOfStates)/2];
        model.setupRateMatrix();
        double[][] Q = model.getRateMatrix();

        // deviations in the detailed balance property (Q_ij*pi_i=Q_ji*pi_j)
        // that applies to time-reversible Q matrices,
        // effectively quantifying the amount of time-nonreversibility in a Q matrix
        int x = 0;
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = i+1; j < nrOfStates; j++) {
                d[x] = Q[i][j]*equilibriumFreqs[i]-Q[j][i]*equilibriumFreqs[j];
                x++;
            }
        }

        return d;
    }

    public void close(PrintStream close) {
    }

    @Override
    public int getDimension() {
        int nrOfStates = this.model.getStateCount();
        return (nrOfStates*nrOfStates - nrOfStates)/2;
    }

    @Override
    public double[] getDoubleValues() {
        return getDeviations();
    }

    @Override
    public double getArrayValue(int dim) {
        return getDoubleValues()[dim];
    }
}
