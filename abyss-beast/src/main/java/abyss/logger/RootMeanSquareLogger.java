package abyss.logger;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.CalculationNode;

import java.io.PrintStream;

import static java.lang.Math.sqrt;

/**
 * @author Jasmine Saghafifar
 */

@Description("Summarises outputs from a logger as a root-mean-square")
public class RootMeanSquareLogger extends CalculationNode implements Loggable, Function {
    final public Input<Function> loggerInput;
    protected Function logger;

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
        double[] values = logger.getDoubleValues();
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
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue(int dim) {
        return getRootMeanSquare();
    }

    @Override
    public double[] getDoubleValues() {
        return null;
    }

}
