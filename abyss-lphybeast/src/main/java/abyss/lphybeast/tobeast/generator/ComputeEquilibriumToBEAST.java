//package abyss.lphybeast.tobeast.generator;
//
//import abyss.ComputeEquilibrium;
//import beast.base.core.BEASTInterface;
//import beast.base.core.Loggable;
//import beast.base.inference.Distribution;
//import beast.base.spec.domain.PositiveReal;
//import beast.base.spec.inference.distribution.Dirichlet;
//import beast.base.spec.inference.distribution.TensorDistribution;
//import beast.base.spec.inference.parameter.RealVectorParam;
//import beast.base.spec.type.Simplex;
//import lphybeast.BEASTContext;
//import lphybeast.GeneratorToBEAST;
//
//import java.util.Arrays;
//
///**
// * @author Jasmine Saghafifar
// */
//public class ComputeEquilibriumToBEAST implements GeneratorToBEAST<ComputeEquilibrium, Distribution> {
//    @Override
//    public Distribution generatorToBEAST(ComputeEquilibrium generator, BEASTInterface value, BEASTContext context) {
//        return createDirichletPrior(generator, value, context);
//    }
//
//    private Distribution createDirichletPrior(ComputeEquilibrium generator, BEASTInterface value, BEASTContext context) {
//        Dirichlet dirichlet = new Dirichlet();
//        double[] a = new double[generator.generate().value().length];
//        Arrays.fill(a, 4.0);
//        RealVectorParam<PositiveReal> alpha = new RealVectorParam<>(a, PositiveReal.INSTANCE);
//        dirichlet.setInputValue("alpha", alpha);
//        dirichlet.initAndValidate();
//        TensorDistribution<Simplex, Double> prior = BEASTContext.createPrior(dirichlet,(RealVectorParam<PositiveReal>) value);
//        context.addExtraLoggable((Loggable) value);
//        return prior;
//    }
//
//    @Override
//    public Class<ComputeEquilibrium> getGeneratorClass() {
//        return ComputeEquilibrium.class;
//    }
//
//    @Override
//    public Class<Distribution> getBEASTClass() {
//        return Distribution.class;
//    }
//}
