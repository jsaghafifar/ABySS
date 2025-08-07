package abyss;

import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Objects;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * @author Jasmine Saghafifar
 */
public class MixedAlignment implements GenerativeDistribution<Alignment> {

    protected static final String aln1ParamName = "aln1";
    protected static final String aln2ParamName = "aln2";
    protected static final String aln1SitesParamName = "sites";
    protected static final String indicatorParamName = "indicator";
    protected static final String siteMixtureWeightsParamName = "weights";

    Value<Alignment> aln1;
    Value<Alignment> aln2;
    Value<Boolean[]> aln1sites;
    Value<Integer> indicator;
    Value<Double[]> weights;

    public MixedAlignment(@ParameterInfo(name = aln1ParamName,
                                  description = "the first simulated alignment.") Value<Alignment> aln1,
                          @ParameterInfo(name = aln2ParamName, description = "the second simulated alignment " +
                                  "(same tree, number of taxa, and length as first).") Value<Alignment> aln2,
                          @ParameterInfo(name = aln1SitesParamName, description = "which sites from the first " +
                                  "alignment will be included.", optional = true) Value<Boolean[]> aln1sites,
                          @ParameterInfo(name = indicatorParamName, description = "the alignment that will be chosen " +
                                  "from aln1 (0), aln2 (1), or a mix (2).") Value<Integer> indicator,
                          @ParameterInfo(name = siteMixtureWeightsParamName, description = "weights that were used to" +
                                  " determine site mixture. Required for lphybeast.",
                                  optional = true) Value<Double[]> weights) {
        this.aln1 = aln1;
        this.aln2 = aln2;

        if (!Objects.equals(aln1.value().nchar(), aln2.value().nchar()))
            LoggerUtils.log.severe("Alignments must have same length.");
        if (aln1.value().ntaxa() != aln2.value().ntaxa())
            LoggerUtils.log.severe("Alignments must have same number of taxa.");

        if (indicator.value() < 3 && indicator.value() > -1)
            this.indicator = indicator;
        else LoggerUtils.log.severe(indicatorParamName+" must be 0, 1, or 2.");

        if (aln1sites.value() != null)
            this.aln1sites = aln1sites;
        else if (indicator.value() == 2)
            LoggerUtils.log.severe("Site mixture model is indicated " +
                    "but no "+aln1SitesParamName+" has been provided.");

        if (weights.value() != null)
            this.weights = weights;
    }

    @GeneratorInfo(name = "MixedAlignment", verbClause = "is created by",
            category = GeneratorCategory.TAXA_ALIGNMENT,
            description = "Simulates model averaging. If indicated, creates a new alignment by mixing two alignments. ")
    @Override
    public RandomVariable<Alignment> sample() {
        // TODO fix viewer bug when indicator comes from a sample array
        final int ind = getModelIndicator().value();
        final Alignment aln1 = getAlignment1().value();
        final Alignment aln2 = getAlignment2().value();

        if (ind==0) return new RandomVariable<>("D", aln1, this);
        else if (ind==1) return new RandomVariable<>("D", aln2, this);
        else {
            int nchar = aln1.nchar();
            int ntaxa = aln1.ntaxa();
            Boolean[] aln1sites = getAlignment1Sites().value();
            Alignment mixAln = new SimpleAlignment(nchar, aln1);
            for (int j = 0; j < nchar; j++) {
                if (aln1sites[j]) {
                    for (int i = 0; i < ntaxa; i++) {
                        int newState = aln1.getState(i, j);
                        mixAln.setState(i, j, newState);
                    }
                } else {
                    for (int i = 0; i < ntaxa; i++) {
                        int newState = aln2.getState(i, j);
                        mixAln.setState(i, j, newState);
                    }
                }
            }
            return new RandomVariable<>("D", mixAln, this);
        }

    }

    public Value<Alignment> getAlignment1() {
        return aln1;
    }

    public Value<Alignment> getAlignment2() {
        return aln2;
    }

    public Value<Boolean[]> getAlignment1Sites() {
        return aln1sites;
    }

    public Value<Integer> getModelIndicator() {
        return indicator;
    }

    public Value<Double[]> getSiteMixtureWeights() {
        return weights;
    }

    @Override
    public SortedMap<String, Value> getParams() {
        SortedMap<String, Value> map = new TreeMap<>();
        map.put(aln1ParamName, aln1);
        map.put(aln2ParamName, aln2);
        if (aln1sites!=null) map.put(aln1SitesParamName, aln1sites);
        map.put(indicatorParamName, indicator);
        if (weights!=null) map.put(siteMixtureWeightsParamName, weights);
        return map;
    }
}
