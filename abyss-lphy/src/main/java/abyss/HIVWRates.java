package abyss;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;

//TODO model citation
public class HIVWRates<U> extends DeterministicFunction<U[]> {

    public HIVWRates() {}

    @GeneratorInfo(name = "hivwRates", description = "")
    public Value<U[]> apply() {
        HIVW hivw = new HIVW(null);
        U[] array = (U[]) hivw.getRates();
        return new Value<>( array, this);
    }
}