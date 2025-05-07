module abyss.lphy {
    requires transitive lphy.base;

    // add any req libs here
//    requires ;

    exports abyss;
    exports abyss.spi;

    // LPhy extensions
    uses lphy.core.spi.Extension;
    // declare what service interface the provider intends to use
    provides lphy.core.spi.Extension with abyss.spi.ABySSImpl;
}