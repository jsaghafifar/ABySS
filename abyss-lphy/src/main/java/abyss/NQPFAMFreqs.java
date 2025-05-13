package abyss;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.Citation;
import lphy.core.model.annotation.GeneratorInfo;

@Citation(value="Dang CC, Minh BQ, McShea H, Masel J, James JE, Vinh LS, Lanfear R. "+
        "nQMaker: Estimating Time Nonreversible Amino Acid Substitution Models, " +
        "Sys Biol, 2022, vol. 71 (pg. 1110-1123)",
        title = "nQMaker: Estimating Time Nonreversible Amino Acid Substitution Models",
        year = 2022,
        authors = {"Dang", "Minh", "McShea", "Masel", "James", "Vinh", "Lanfear"},
        DOI = "https://doi.org/10.1093/sysbio/syac007")
public class NQPFAMFreqs<U> extends DeterministicFunction<U[]> {

    public NQPFAMFreqs() {}

    @GeneratorInfo(name = "nqpfamRates", description = "")
    public Value<U[]> apply() {
        NQPFAM nqpfam = new NQPFAM(null);
        U[] array = (U[]) nqpfam.getFreq();
        return new Value<>( array, this);
    }

}