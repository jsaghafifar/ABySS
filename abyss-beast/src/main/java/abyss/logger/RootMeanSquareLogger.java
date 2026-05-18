package abyss.logger;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.RealVector;

import java.io.PrintStream;

import static java.lang.Math.sqrt;

/**
 * @author Jasmine Saghafifar
 */

@Description("Summarises outputs from a logger as a root-mean-square")
public class RootMeanSquareLogger extends CalculationNode implements Loggable, RealScalar<NonNegativeReal> {
    final public Input<RealVector<?>> loggerInput;
    protected RealVector<?> logger;

    public RootMeanSquareLogger() {
        this.loggerInput = new Input<>("logger", "Logger to be summarised.",
                Input.Validate.REQUIRED);
    }

    public void initAndValidate() {
        this.logger = this.loggerInput.get();
    }

    public void init(PrintStream out) {
        out.print(getID() + "\t");
    }

    public void log(long sample, PrintStream out) {
        out.print(getRootMeanSquare() + "\t");
    }

    private double getRootMeanSquare() {
        Double[] values = logger.getElements().toArray(new Double[0]);
        double rms = 0;
        for (double value : values) {
            rms += value * value;
        }
        rms /= values.length;
        rms = sqrt(rms);

        return rms;
    }

    public void close(PrintStream close) {
    }

    @Override
    public int size() {
        return 1;
    }

    @Override
    public double get() {
        return getRootMeanSquare();
    }

    @Override
    public NonNegativeReal getDomain() {
        return NonNegativeReal.INSTANCE;
    }
}
