set(LIBINTINTEGRALS_HEADERS
        LibintIntegrals/LibintIntegrals.h
        LibintIntegrals/Libint.h
        LibintIntegrals/BasisSetHandler.h
        LibintIntegrals/OneBodyIntegrals.h
        LibintIntegrals/IntegralEvaluatorSettings.h
        LibintIntegrals/TwoBodyIntegrals/Digester.h
        LibintIntegrals/TwoBodyIntegrals/Evaluator.h
        LibintIntegrals/TwoBodyIntegrals/Prescreener.h
        LibintIntegrals/TwoBodyIntegrals/SaverDigester.h
        LibintIntegrals/TwoBodyIntegrals/SymmetryHelper.h
        LibintIntegrals/TwoBodyIntegrals/VoidPrescreener.h
        LibintIntegrals/TwoBodyIntegrals/COMSaverDigester.h
        LibintIntegrals/TwoBodyIntegrals/CauchySchwarzDensityPrescreener.h
        LibintIntegrals/TwoBodyIntegrals/HartreeFock/CoulombExchangeDigester.h
        LibintIntegrals/TwoBodyIntegrals/HartreeFock/CoulombExchangeConstructor.h
        LibintIntegrals/TwoBodyIntegrals/HartreeFock/TwoTypeCoulombDigester.h
        LibintIntegrals/TwoBodyIntegrals/HartreeFock/TwoTypeCoulombConstructor.h
        )

set(LIBINTINTEGRALS_SOURCES
        LibintIntegrals/LibintIntegrals.cpp
        LibintIntegrals/BasisSetHandler.cpp
        LibintIntegrals/OneBodyIntegrals.cpp
        LibintIntegrals/Libint.cpp
        LibintIntegrals/TwoBodyIntegrals/SaverDigester.cpp
        LibintIntegrals/TwoBodyIntegrals/COMSaverDigester.cpp
        LibintIntegrals/TwoBodyIntegrals/CauchySchwarzDensityPrescreener.cpp
        LibintIntegrals/TwoBodyIntegrals/HartreeFock/CoulombExchangeDigester.cpp
        LibintIntegrals/TwoBodyIntegrals/HartreeFock/CoulombExchangeConstructor.cpp
        LibintIntegrals/TwoBodyIntegrals/HartreeFock/TwoTypeCoulombDigester.cpp
        LibintIntegrals/TwoBodyIntegrals/HartreeFock/TwoTypeCoulombConstructor.cpp
        )