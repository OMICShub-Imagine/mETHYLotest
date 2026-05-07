.onAttach <- function(libname, pkgname) {
  v <- utils::packageVersion(pkgname)
  packageStartupMessage(
    "\n",
    "            ______ _______ _    ___     ___           _            _   \n",
    "           |  ____|__   __| |  | \\ \\   / / |         | |          | |  \n",
    "  _ __ ___ | |__     | |  | |__| |\\ \\_/ /| |     ___ | |_ ___  ___| |_ \n",
    " | '_ ` _ \\|  __|    | |  |  __  | \\   / | |    / _ \\| __/ _ \\/ __| __|\n",
    " | | | | | | |____   | |  | |  | |  | |  | |___| (_) | ||  __/\\__ \\ |_ \n",
    " |_| |_| |_|______|  |_|  |_|  |_|  |_|  |______\\___/ \\__\\___||___/\\__|\n",
    "\n",
    "  ===========================================================\n",
    "    mETHYLotest v", v, "\n",
    "    DNA Methylation Analysis Pipeline\n",
    "  ===========================================================\n",
    "\n",
    "    EPIC module   -   Illumina 450K / EPICv1 / EPICv2 / Mouse\n",
    "    NGS module   -   WGBS / RRBS / Nanopore / PacBio\n",
    "\n",
    "  ===========================================================\n",
    "    mETHYLotest.EPIC.pipeline()   EPIC array analysis\n",
    "    mETHYLotest.NGS.pipeline()    NGS methylation analysis\n",
    "  ===========================================================\n"
  )
}
