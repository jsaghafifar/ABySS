package abyss.logger;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.type.Simplex;
import beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.DoubleStream;

/**
 * @author Jasmine Saghafifar
 */

@Description("Equilibrium frequency logger")
public class ClassicFrequencyLogger extends CalculationNode implements Loggable, Simplex {
    final public Input<SVSGeneralSubstitutionModel> modelInput;
    final public Input<String> keysInput;
    protected SVSGeneralSubstitutionModel model;
    protected List<String> keys;
    private static final double DEFAULT_BRANCH_LENGTH = 100000;

    public ClassicFrequencyLogger() {
        this.modelInput = new Input<>("model", "ABYSS SVS general substitution model.", Input.Validate.REQUIRED);
        this.keysInput = new Input<>("keys", "freq state names");
    }

    public void initAndValidate() {
        this.model = this.modelInput.get();

        if (this.keysInput.get() != null) {
            String[] keysArr = keysInput.get().split(" ");
            List<String> keys = Collections.unmodifiableList(Arrays.asList(keysArr));
            if (keys.size() != size())
                throw new IllegalArgumentException("For 1D array, keys must have the same length as dimension ! " +
                        "Dimension = " + size() + ", but keys.size() = " + keys.size());
            this.keys = keys;
        }
    }

    public void init(PrintStream out) {
        // substModel + freq + state name of j = param name
        String id = getID();
        String stateName;
        for (int i = 0; i < this.model.getStateCount(); i++) {
            stateName = this.keys != null ? keys.get(i) : String.valueOf(i);
            out.print(id + '.' + stateName + "\t");
        }

    }

    public void log(long sample, PrintStream out) {
        for (int i = 0; i < this.model.getStateCount(); i++) {
            out.print(get(i) + "\t");
        }

    }

    private List<Double> getFrequencies(SVSGeneralSubstitutionModel model) {
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
        return DoubleStream.of(f).boxed().toList();
    }

    public void close(PrintStream close) {
    }

    @Override
    public int size() {
        return this.model.getStateCount();
    }

    @Override
    public List<Double> getElements() {
        return getFrequencies(this.model);
    }

    @Override
    public double get(int i) {
        return getElements().get(i);
    }

    @Override
    public UnitInterval getDomain() {
        return UnitInterval.INSTANCE;
    }
}
