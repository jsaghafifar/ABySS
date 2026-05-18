package abyss.substitutionmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.substitutionmodel.ComplexColtEigenSystem;
import beast.base.spec.domain.NonNegativeInt;
import beast.base.spec.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.spec.inference.parameter.IntScalarParam;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Jasmine Saghafifar
 */

@Description("Substitution model that can average over a number of substitution models ")

public class ABySSModelAveraging extends GeneralSubstitutionModel {
    final public Input<List<GeneralSubstitutionModel>> substModelInput =
            new Input<>("model", "substitution model", new ArrayList<>(), Input.Validate.REQUIRED);
    final public Input<IntScalarParam<NonNegativeInt>> modelIndicatorInput =
            new Input<>("modelIndicator", "index of the model in list of models that is used for its rates and frequencies", Input.Validate.REQUIRED);
    final public Input<Integer> nrOfStatesInput =
            new Input<>("nrOfStates", "number of states parameter", Input.Validate.REQUIRED);

    IntScalarParam<NonNegativeInt> modelIndicator;
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

        updateMatrix = true;
        eigenSystem = new ComplexColtEigenSystem(nrOfStates);

        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[nrOfStates*(nrOfStates-1)];
        storedRelativeRates = new double[nrOfStates*(nrOfStates-1)];
    }

    @Override
    public void setupRelativeRates() {
        GeneralSubstitutionModel model = models.get(modelIndicator.get());
        model.setupRelativeRates();
        relativeRates = model.getRelativeRates();
    }

    @Override
    public void setupRateMatrix() {
        GeneralSubstitutionModel model = models.get(modelIndicator.get());
        model.setupRateMatrix();
        rateMatrix = model.getRateMatrix();
    }

    @Override
    public double[] getFrequencies() {
        GeneralSubstitutionModel model = models.get(modelIndicator.get());
        return model.getFrequencies();
    }

    @Override
    public double[][] getRateMatrix() {
        GeneralSubstitutionModel model = models.get(modelIndicator.get());
        model.setupRateMatrix();
        rateMatrix = model.getRateMatrix();
        return rateMatrix;
    }

    @Override
    public double[] getRateMatrix(Node node) {
        GeneralSubstitutionModel model = models.get(modelIndicator.get());
        model.setupRateMatrix();
        double[][] matrix = model.getRateMatrix();
        double[] rates = new double[nrOfStates*nrOfStates];
        for (int i = 0; i < nrOfStates; i++) {
            // TODO check array copy done right. was i*nrOfStates + j = i,j
            System.arraycopy(matrix[i], 0, rates, i * nrOfStates, nrOfStates);
        }
        return rates;
    }

    @Override
    protected boolean requiresRecalculation() {
        if (modelIndicator.somethingIsDirty()) {
            updateMatrix = true;
            return true;
        }
        return super.requiresRecalculation();
    }
}
