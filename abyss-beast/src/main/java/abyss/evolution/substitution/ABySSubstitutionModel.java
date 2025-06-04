package abyss.evolution.substitution;


import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.substitutionmodel.ComplexSubstitutionModel;
import beast.base.evolution.substitutionmodel.DefaultEigenSystem;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.Parameter;
import abyss.inference.AbyssSVS;

import java.util.Arrays;

//TODO Add checks for connected graph?

/**
 * ported from beast1 to BEAST_CLASSIC - author: Marc Suchard
 * @author dkuh004
 *         Date: Sep 18, 2011
 *         Time: 6:03:49 PM
 */
@Description("ABySS ver SVS General Substitution Model")
public class ABySSubstitutionModel extends ComplexSubstitutionModel implements AbyssSVS {

    public Input<BooleanParameter> indicator = new Input<BooleanParameter>("rateIndicator",
            "rates to indicate the presence or absence of transition matrix entries", Validate.REQUIRED);

    public Input<Boolean> isSymmetricInput = new Input<Boolean>("symmetric",
    		"Indicates the rate matrix is symmetric. " +
            "If true (default) n(n-1)/2 rates and indicators need to be specified. " +
            "If false, n(n-1) rates and indicators need to be specified.", Boolean.TRUE);

    private BooleanParameter rateIndicator;
    private boolean isSymmetric = false;
    private static final double BRANCH_LENGTH = 1e6;

    @Override
    public void initAndValidate(){

        frequencies = frequenciesInput.get();

        updateMatrix = true;
        nrOfStates = frequencies.getFreqs().length;
        if (isSymmetricInput.get() && ratesInput.get().getDimension() != nrOfStates * (nrOfStates-1)/2) {
            throw new IllegalArgumentException("Dimension of input 'rates' is " + ratesInput.get().getDimension() + " but a " +
                    "rate matrix of dimension " + nrOfStates + "x" + (nrOfStates -1) + "/2" + "=" + nrOfStates * (nrOfStates -1) / 2 + " was " +
                    "expected");
        }
        if (!isSymmetricInput.get() && ratesInput.get().getDimension() != nrOfStates * (nrOfStates-1)) {
            int dim = nrOfStates * (nrOfStates -1);
            Log.warning.println("Dimension of input 'rates' is " + ratesInput.get().getDimension() + ". " +
            		"Changing dimension to " + nrOfStates + "x" + (nrOfStates -1)  + "=" + dim);
            if (ratesInput.get() instanceof Parameter.Base) {
            	((Parameter.Base)ratesInput.get()).setDimension(dim);
            }
            Log.warning.println("Setting dimension of indicators to " + dim);
            indicator.get().setDimension(dim);
            isSymmetric = false;
        }

        if (!isSymmetricInput.get() && eigenSystemClass.equals(DefaultEigenSystem.class.getName())) {
        	Log.warning.println("WARNING: eigenSystemClass is DefautlEigneSystem, which may cause trouble with asymtric analysis. "
        			+ "You may want to consider eigensystem='beast.evolution.substitutionmodel.RobustEigenSystem' instead.");
        }

        eigenSystem = createEigenSystem();
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[ratesInput.get().getDimension()];
        storedRelativeRates = new double[ratesInput.get().getDimension()];

        rateIndicator = indicator.get();
        super.initAndValidate();
    }


    public Parameter<?> getIndicators() {
        return rateIndicator;
    }

    public boolean validState() {
        return !updateMatrix || Utils.connectedAndWellConditioned(probability,this);
    }


    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        updateMatrix = true;
    }

    private double[] probability = null;

    @Override
    public void setupRelativeRates() {

        Function rates = this.ratesInput.get();
        for (int i = 0; i < relativeRates.length; i++) {
            relativeRates[i] = rates.getArrayValue(i) * (rateIndicator.getValue(i)?1.:0.);
        }
    }

    /** sets up rate matrix **/
    @Override
    public void setupRateMatrix() {
        double [] fFreqs = frequencies.getFreqs();
        if (isSymmetricInput.get()) {
            int count = 0;
            for (int i = 0; i < nrOfStates; i++) {
                rateMatrix[i][i] = 0;
                for (int j = i+1; j < nrOfStates; j++) {
                    rateMatrix[i][j] = relativeRates[count];
                    rateMatrix[j][i] = relativeRates[count];
                    count++;
                }
            }
            // bring in frequencies
            for (int i = 0; i < nrOfStates; i++) {
                for (int j = i + 1; j < nrOfStates; j++) {
                    rateMatrix[i][j] *= fFreqs[j];
                    rateMatrix[j][i] *= fFreqs[i];
                }
            }
    	} else {
            for (int i = 0; i < nrOfStates; i++) {
                rateMatrix[i][i] = 0;
                for (int j = 0; j < i; j++) {
                    rateMatrix[i][j] = relativeRates[i * (nrOfStates - 1) + j];
                }
                for (int j = i + 1; j < nrOfStates; j++) {
                    rateMatrix[i][j] = relativeRates[i * (nrOfStates - 1) + j - 1];
                }
            }
        }

        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double fSum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    fSum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -fSum;
        }

        // normalise rate matrix to one expected substitution per unit time
        double[] f;
        if (!isSymmetric) {
            // compute equilibrium frequencies if nonreversible
            f = getEquilibriumFrequencies(rateMatrix);
        } else f = fFreqs;
        double fSubst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            fSubst += -rateMatrix[i][i] * f[i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / fSubst;
            }
        }
    } // setupRateMatrix


    @Override
    protected boolean requiresRecalculation() {
    	return super.requiresRecalculation();
    }
    
    @Override
    public boolean canReturnComplexDiagonalization() {
    	return !isSymmetric;
    }

    @Override
    public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
        try {super.getTransitionProbabilities( node,  startTime,  endTime,  rate, matrix);} catch (Exception exception) {
            Arrays.fill(matrix, 0.0);
        }

    }

    private double[] getEquilibriumFrequencies(double[][] Qm) {
        double temp;

        EigenDecomposition decomposition = eigenSystem.decomposeMatrix(Qm);
        double[] Eval = decomposition.getEigenValues();
        double[] evec = decomposition.getEigenVectors();
        double[] ievc = decomposition.getInverseEigenVectors();

        double[][] iexp = new double[nrOfStates][nrOfStates];
        double[] EvalImag = new double[nrOfStates];
        double[][] transProbs = new double[nrOfStates][nrOfStates];

        System.arraycopy(Eval, nrOfStates, EvalImag, 0, nrOfStates);
        for (int i = 0; i < nrOfStates; i++) {
            if (EvalImag[i] == 0) {
                // 1x1 block
                temp = Math.exp(BRANCH_LENGTH * Eval[i]);
                for (int j = 0; j < nrOfStates; j++) {
                    iexp[i][j] = ievc[i * nrOfStates + j] * temp;
                }
            } else {
                // 2x2 conjugate block
                // If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
                // exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
                int i2 = i + 1;
                double b = EvalImag[i];
                double expat = Math.exp(BRANCH_LENGTH * Eval[i]);
                double expatcosbt = expat * Math.cos(BRANCH_LENGTH * b);
                double expatsinbt = expat * Math.sin(BRANCH_LENGTH * b);

                for (int j = 0; j < nrOfStates; j++) {
                    iexp[i][j] = expatcosbt * ievc[i * nrOfStates + j] +
                            expatsinbt * ievc[i2 * nrOfStates + j];
                    iexp[i2][j] = expatcosbt * ievc[i2 * nrOfStates + j] -
                            expatsinbt * ievc[i * nrOfStates + j];
                }
                i++; // processed two conjugate rows
            }
        }

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                temp = 0.0;
                for (int k = 0; k < nrOfStates; k++) {
                    temp += evec[i * nrOfStates + k] * iexp[k][j];
                }
                transProbs[i][j] = Math.abs(temp);
            }
        }
        double[] freqs = new double[transProbs.length];
        for (int i = 0; i < freqs.length; i++) {
            freqs[i] = transProbs[0][i];
            for (int j = 1; j < freqs.length; j++) {
//                if (Math.abs(transProbs[0][i] - transProbs[j][i]) > 1e-6) {
//                    System.out.println("WARNING: branch length used to get equilibrium distribution was not long enough!");
//                    break;
//                }

            }
        }
        return freqs;
    }

}
