module abyss.lphy {
    requires transitive lphy.base;
    requires colt;

    exports abyss.lphy;

    // LPhy extensions
    uses lphy.core.spi.Extension;
    // declare what service interface the provider intends to use
    provides lphy.core.spi.Extension with abyss.spi.ABySSImpl;
}