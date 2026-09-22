open module abyss.lphybeast {
    requires lphy.beast;
    requires abyss.lphy;
    requires abyss.beast;
    requires beast.base;
    requires lphy.base;
    requires beast.classic;
    requires beast.pkgmgmt;

    exports abyss.lphybeast.spi;
    exports abyss.lphybeast.tobeast.generator;

    provides lphybeast.spi.LPhyBEASTMapping with abyss.lphybeast.spi.ABySSLBExtImpl;
}