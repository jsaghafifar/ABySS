package abyss.evolution.substitution;


import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.substitutionmodel.ComplexSubstitutionModel;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.Parameter;
import abyss.inference.AbyssSVS;

import java.util.Arrays;


/**
 * ported from beast1 to BEAST_CLASSIC - author: Marc Suchard
 * @author dkuh004
 *         Date: Sep 18, 2011
 *         Time: 6:03:49 PM
 */
@Description("ABySS ver SVS General Substitution Model")
public class ABySSubstitutionModel extends ComplexSubstitutionModel implements AbyssSVS {

    public Input<BooleanParameter> indicator = new Input<BooleanParameter>("rateIndicator",
            "rates to indicate the presence or absence of transition matrix entries", Validate.OPTIONAL);

    public Input<Boolean> isSymmetricInput = new Input<Boolean>("symmetric",
    		"Indicates the rate matrix is symmetric. " +
            "If true (default) n(n-1)/2 rates and indicators need to be specified. " +
            "If false, n(n-1) rates and indicators need to be specified.", Boolean.FALSE);

    public ABySSubstitutionModel() {
        frequenciesInput.setRule(Validate.OPTIONAL);
    }

    private BooleanParameter rateIndicator;
    private boolean isSymmetric = false;
    private static final double DEFAULT_BRANCH_LENGTH = 100000;

    @Override
    public void initAndValidate(){
        if (indicator.get() != null) {
            rateIndicator = indicator.get();

            if (indicator.get().getDimension() != ratesInput.get().getDimension())
                throw new IllegalArgumentException("Dimension of inputs 'rates' and 'rateIndicator' must match.");
        }

        if (isSymmetricInput.get()) {
            frequencies = frequenciesInput.get();
            nrOfStates = frequencies.getFreqs().length;

            if (ratesInput.get().getDimension() != nrOfStates * (nrOfStates-1)/2) {
                throw new IllegalArgumentException("Dimension of input 'rates' is " + ratesInput.get().getDimension() + " but a " +
                        "rate matrix of dimension " + nrOfStates + "x" + (nrOfStates -1) + "/2" + "=" + nrOfStates * (nrOfStates -1) / 2 + " was " +
                        "expected");
            }
        } else {
            if (frequenciesInput.get() != null)
                throw new IllegalArgumentException("Frequencies cannot be specified for nonreversible substitution");
            double root = (-1 - Math.sqrt(1+4* ratesInput.get().getDimension()))/2;
            if (root - (int) root < 1e-6) nrOfStates = (int) Math.abs(root);
            else throw new IllegalArgumentException("Rates must have dimension nrOfStates*(nrOfStates-1) when nonreversible");

        }

        updateMatrix = true;
        eigenSystem = createEigenSystem(); // ComplexColtEigenSystem
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[ratesInput.get().getDimension()];
        storedRelativeRates = new double[ratesInput.get().getDimension()];

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
        if (this.indicator.get() != null) {
            for (int i = 0; i < relativeRates.length; i++) {
                relativeRates[i] = rates.getArrayValue(i) * (rateIndicator.getValue(i)?1.:0.);
            }
        } else for (int i = 0; i < relativeRates.length; i++) {
            relativeRates[i] = rates.getArrayValue(i);
        }
    }

    /** sets up rate matrix **/
    @Override
    public void setupRateMatrix() {
        double[] f;
        if (isSymmetricInput.get()) {
            f = frequencies.getFreqs();
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
                    rateMatrix[i][j] *= f[j];
                    rateMatrix[j][i] *= f[i];
                }
            }
    	} else { // nonreversible Q mat not constructed from freqs
            f = new double[nrOfStates];
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
        if (!isSymmetricInput.get()) f = this.getFrequencies();

        // normalise rate matrix to one expected substitution per unit time
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

    @Override
    public double[] getFrequencies() {
        if (!isSymmetric) {
            setupRelativeRates();
            double[][] Qm = new double[nrOfStates][nrOfStates];
            for (int i = 0; i < nrOfStates; i++) {
                Qm[i][i] = 0;
                for (int j = 0; j < i; j++) {
                    Qm[i][j] = relativeRates[i * (nrOfStates - 1) + j];
                }
                for (int j = i + 1; j < nrOfStates; j++) {
                    Qm[i][j] = relativeRates[i * (nrOfStates - 1) + j - 1];
                }
            }
            for (int i = 0; i < nrOfStates; i++) {
                double fSum = 0.0;
                for (int j = 0; j < nrOfStates; j++) {
                    if (i != j)
                        fSum += Qm[i][j];
                }
                Qm[i][i] = -fSum;
            }
            return getEquilibriumFrequencies(Qm);
        }
        return this.frequencies.getFreqs();
    }

    private double[] getEquilibriumFrequencies(double[][] Qm) {
        double temp;

        EigenDecomposition complexDecomposition = eigenSystem.decomposeMatrix(Qm);
        double[] Eval = complexDecomposition.getEigenValues();
        double[] evec = complexDecomposition.getEigenVectors();
        double[] ievc = complexDecomposition.getInverseEigenVectors();

        double[][] iexp = new double[nrOfStates][nrOfStates];
        double[] EvalImag = new double[nrOfStates];
        double[][] transProbs = new double[nrOfStates][nrOfStates];
        double[] freqs = new double[transProbs.length];

        System.arraycopy(Eval, nrOfStates, EvalImag, 0, nrOfStates);
        double t = DEFAULT_BRANCH_LENGTH;
        boolean equilibrium = false;

        while (!equilibrium) {
            for (int i = 0; i < nrOfStates; i++) {
                if (EvalImag[i] == 0) {
                    // 1x1 block
                    temp = Math.exp(t * Eval[i]);
                    for (int j = 0; j < nrOfStates; j++) {
                        iexp[i][j] = ievc[i * nrOfStates + j] * temp;
                    }
                } else {
                    // 2x2 conjugate block
                    // If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
                    // exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
                    int i2 = i + 1;
                    double b = EvalImag[i];
                    double expat = Math.exp(t * Eval[i]);
                    double expatcosbt = expat * Math.cos(t * b);
                    double expatsinbt = expat * Math.sin(t * b);

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

            boolean reached = true;
            for (int i = 0; i < nrOfStates; i++) {
                freqs[i] = transProbs[0][i];
                for (int j = 1; j < nrOfStates; j++) {
                    if (Math.abs(transProbs[0][i] - transProbs[j][i]) > 1e-6) {
                        reached = false;
                        break;
                    }
                }
                if (!reached) break;
            }
            if (reached) equilibrium = true;
            t *= 10;
        }
        return freqs;
    }

}
