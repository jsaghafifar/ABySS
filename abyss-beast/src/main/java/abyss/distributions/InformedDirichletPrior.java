package abyss.distributions;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;


@Description("Informed Dirichlet distribution.")
public class InformedDirichletPrior extends ParametricDistribution {
    final public Input<Function> alphaInput = new Input<>("alpha", "coefficients of the Dirichlet distribution", Validate.REQUIRED);
    final public Input<RealParameter> scaleInput = new Input<>("scale", "scale of coefficients", Validate.REQUIRED);

    @Override
    public void initAndValidate() {
    }

    @Override
    public Distribution getDistribution() {
        return null;
    }

    class DirichletImpl implements ContinuousDistribution {
        Double[] m_fAlpha;

        void setAlpha(Double[] alpha) {
            m_fAlpha = alpha;
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
            throw new MathException("Not implemented yet");
        }

        @Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            throw new MathException("Not implemented yet");
        }

        @Override
        public double inverseCumulativeProbability(double p) throws MathException {
            throw new MathException("Not implemented yet");
        }

        @Override
        public double density(double x) {
            return Double.NaN;
        }

        @Override
        public double logDensity(double x) {
            return Double.NaN;
        }
        
    } // class DirichletImpl


    @Override
    public double calcLogP(Function pX) { //TODO add scale and fixed minimum here
        double[] alpha = alphaInput.get().getDoubleValues();
        if (alphaInput.get().getDimension() != pX.getDimension()) {
            throw new IllegalArgumentException("Dimensions of alpha and x should be the same, but dim(alpha)=" + alphaInput.get().getDimension()
                    + " and dim(x)=" + pX.getDimension());
        }
        double logP = 0;
        double sumAlpha = 0;
        for (int i = 0; i < pX.getDimension(); i++) {
            double x = pX.getArrayValue(i);
            logP += (alpha[i] - 1) * Math.log(x);
            logP -= org.apache.commons.math.special.Gamma.logGamma(alpha[i]);
            sumAlpha += alpha[i];
        }
        logP += org.apache.commons.math.special.Gamma.logGamma(sumAlpha);
        return logP;
    }

	@Override
	public Double[][] sample(int size) {
		int dim = alphaInput.get().getDimension();

        Double[] conc = new Double[dim];
        double concSum = 0.0;
        for (int i = 0; i < dim; i++) {
            conc[i] = alphaInput.get().getArrayValue(i);
            concSum += conc[i];
        }

        Double[][] samples = new Double[size][];
        double scale = scaleInput.get().getValue();
        double min = scale * dim / 200000;
		for (int i = 0; i < size; i++) {
			Double[] dirichletSample = new Double[dim];
			double dirSum = 0.0;
            int x = 0;
			for (int j = 0; j < dim; j++) {
                conc[j] /= concSum;
                conc[j] *= scale * dim;
				dirichletSample[j] = Randomizer.nextGamma(conc[j], 1.0);
                if (dirichletSample[j] < min) {dirichletSample[j] = min; x++;} //TODO add check for >25% samples fixed to minimum
				dirSum += dirichletSample[j];
			}
			for (int j = 0; j < dim; j++) {
				dirichletSample[j] = dirichletSample[j] / dirSum;
			}
			samples[i] = dirichletSample;

		}
		return samples;
	}
}
