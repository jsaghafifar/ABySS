package abyss.lphybeast.tobeast.generator;

import abyss.ComputeEquilibrium;
import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.core.Loggable;
import beast.base.inference.distribution.Dirichlet;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.RealParameter;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.Arrays;

public class ComputeEquilibriumToBEAST implements GeneratorToBEAST<ComputeEquilibrium, Prior> {
    @Override
    public Prior generatorToBEAST(ComputeEquilibrium generator, BEASTInterface value, BEASTContext context) {
        return createDirichletPrior(generator, value, context);
    }

    private Prior createDirichletPrior(ComputeEquilibrium generator, BEASTInterface value, BEASTContext context) {
        Dirichlet dirichlet = new Dirichlet();
        Double[] a = new Double[generator.generate().value().length];
        Arrays.fill(a, 4.0);
        Function alpha = new RealParameter(a);
        dirichlet.setInputValue("alpha", alpha);
        dirichlet.initAndValidate();
        Prior prior = BEASTContext.createPrior(dirichlet,(RealParameter) value);
        context.addExtraLoggable((Loggable) value);
        return prior;
    }

    @Override
    public Class<ComputeEquilibrium> getGeneratorClass() {
        return ComputeEquilibrium.class;
    }

    @Override
    public Class<Prior> getBEASTClass() {
        return Prior.class;
    }
}
