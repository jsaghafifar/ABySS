package abyss.logger;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.Real;
import beast.base.spec.type.RealVector;
import beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.DoubleStream;

/**
 * @author Jasmine Saghafifar
 */

@Description("Logger for deviations from detailed balance in nonreversible Q")
public class ClassicDetailedBalanceLogger extends CalculationNode implements Loggable, RealVector<Real> {
    final public Input<SVSGeneralSubstitutionModel> modelInput;
    final public Input<String> keysInput;
    protected SVSGeneralSubstitutionModel model;
    protected List<String> keys;
    protected double[] equilibriumFreqs;

    public ClassicDetailedBalanceLogger() {
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
            if (keys.size() != size())
                throw new IllegalArgumentException("Number of keys must match dimension. " +
                        "Dimension = " + size() + ", but keys.size() = " + keys.size());
            this.keys = keys;
        }
    }

    public void init(PrintStream out) {
        // logger id + rate name = param name
        String id = getID();
        String rateName;
        for (int i = 0; i < size(); i++) {
            rateName = this.keys != null ? keys.get(i) : String.valueOf(i);
            out.print(id + '.' + rateName + "\t");
        }

    }

    public void log(long sample, PrintStream out) {
        for (int i = 0; i < size(); i++) {
            out.print(get(i) + "\t");
        }

    }

    private List<Double> getDeviations() {
        int nrOfStates = model.getStateCount();
        double[] d = new double[(nrOfStates*nrOfStates-nrOfStates)/2];
        model.setupRateMatrix();
        this.equilibriumFreqs = model.getFrequencies();
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

        return DoubleStream.of(d).boxed().toList();
    }

    public void close(PrintStream close) {
    }

    @Override
    public int size() {
        int nrOfStates = this.model.getStateCount();
        return (nrOfStates*nrOfStates - nrOfStates)/2;
    }

    @Override
    public double get(int i) {
        return getElements().get(i);
    }

    @Override
    public List<Double> getElements() {
        return getDeviations();
    }

    @Override
    public Real getDomain() {
        return Real.INSTANCE;
    }
}
