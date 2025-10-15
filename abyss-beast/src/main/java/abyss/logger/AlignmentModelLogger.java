package abyss.logger;

import abyss.distributions.MixedTreeLikelihood;
import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.util.Randomizer;
import org.apache.commons.math3.util.FastMath;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * @author Jasmine Saghafifar
 */
@Description("Model logger for models with a number of mixed likelihoods.")
public class AlignmentModelLogger extends BEASTObject implements Loggable {

    final public Input<MixedTreeLikelihood> mixedLikelihoodsInput = new Input<>("mixedLikelihoods",
            "Mixed tree likelihoods", Input.Validate.REQUIRED);

    int modelCount;

    @Override
    public void initAndValidate() {
        this.modelCount = mixedLikelihoodsInput.get().pLikelihoods.get().size();
    }


    @Override
    public void init(PrintStream out) {
        out.print((getID() != null ? getID() : "modelIndicator") + "\t");
    }

    @Override
    public void log(long sample, PrintStream out) {
        // sample from models posteriors
        double[] partialLogLikelihoods = mixedLikelihoodsInput.get().getPartialLogLikelihoods();
        double[] posteriorOfEachModel = new double[partialLogLikelihoods.length];

        double max = Arrays.stream(partialLogLikelihoods).max().getAsDouble();
        double pSum = 0;
        for (int i = 0; i < partialLogLikelihoods.length; i++) {
            posteriorOfEachModel[i] = FastMath.exp(partialLogLikelihoods[i]-max);
            pSum += posteriorOfEachModel[i];
        }
        for (int i = 0; i < posteriorOfEachModel.length; i++) {
            posteriorOfEachModel[i] /= pSum;
        }

        int modelIndex = Randomizer.randomChoicePDF(posteriorOfEachModel);
        out.print(modelIndex + "\t");

    }

    @Override
    public void close(PrintStream out) {

    }
}
