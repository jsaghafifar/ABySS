package abyss.substitutionmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.substitutionmodel.ComplexColtEigenSystem;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.IntegerParameter;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.List;

@Description("Substitution model that can average over a number of substitution models ")

public class ABySSModelAveraging extends GeneralSubstitutionModel {
    final public Input<List<GeneralSubstitutionModel>> substModelInput =
            new Input<>("model", "substitution model", new ArrayList<>(), Input.Validate.REQUIRED);
    final public Input<IntegerParameter> modelIndicatorInput =
            new Input<>("modelIndicator", "index of the model in list of models that is used for its rates and frequencies", Input.Validate.REQUIRED);
    final public Input<Integer> nrOfStatesInput =
            new Input<>("nrOfStates", "number of states parameter", Input.Validate.REQUIRED);

    IntegerParameter modelIndicator;
    List<GeneralSubstitutionModel> models;

    public ABySSModelAveraging() {
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
        ratesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        nrOfStates = nrOfStatesInput.get();
        if (frequenciesInput.get() != null) {
            frequencies = frequenciesInput.get();
            if (frequencies.getFreqs().length != nrOfStates) throw new IllegalArgumentException("Wrong nrOfStates or freq dim");
        }

        models = substModelInput.get();
        modelIndicator = modelIndicatorInput.get();
        if (modelIndicator.getUpper() > models.size() - 1) {
            Log.warning("Setting upper limit of " + modelIndicator.getID() + " to " + (models.size()-1) +".");
            modelIndicator.setUpper(models.size() - 1);
        }
        if (modelIndicator.getLower() < 0) {
            Log.warning("Setting lower limit of " + modelIndicator.getID() + " to 0.");
            modelIndicator.setLower(0);
        }

        updateMatrix = true;
        eigenSystem = new ComplexColtEigenSystem(nrOfStates);

        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[nrOfStates*(nrOfStates-1)];
        storedRelativeRates = new double[nrOfStates*(nrOfStates-1)];
    }

    @Override
    public void setupRelativeRates() {
        GeneralSubstitutionModel model = models.get(modelIndicator.getValue());
        model.setupRelativeRates();
        double[] rates = model.getRelativeRates();
        relativeRates = rates;
    }

    @Override
    public void setupRateMatrix() {
        GeneralSubstitutionModel model = models.get(modelIndicator.getValue());
        model.setupRateMatrix();
        double[][] matrix = model.getRateMatrix();
        rateMatrix = matrix;
    }

    @Override
    public double[] getFrequencies() {
        GeneralSubstitutionModel model = models.get(modelIndicator.getValue());
        double[] f = model.getFrequencies();
        return f;
    }

    @Override
    public double[][] getRateMatrix() {
        GeneralSubstitutionModel model = models.get(modelIndicator.getValue());
        model.setupRateMatrix();
        double[][] matrix = model.getRateMatrix();
        rateMatrix = matrix;
        return rateMatrix;
    }

    @Override
    public double[] getRateMatrix(Node node) {
        GeneralSubstitutionModel model = models.get(modelIndicator.getValue());
        model.setupRateMatrix();
        double[][] matrix = model.getRateMatrix();
        double[] rates = new double[nrOfStates*nrOfStates];
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rates[i*nrOfStates + j] = matrix[i][j];
            }
        }
        return rates;
    }

    @Override
    protected boolean requiresRecalculation() {
        if (modelIndicator.isDirtyCalculation()) {
            updateMatrix = true;
            return true;
        }
        return super.requiresRecalculation();
    }
}
