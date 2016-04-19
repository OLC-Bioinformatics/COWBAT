#!/usr/bin/env python
from Bio.Application import _Option, AbstractCommandline, _Switch, _Argument
import re

__author__ = 'mike knowles'
__doc__ = 'Wrapper for bowtie2'


class _PipeArgumentList(_Argument):
    """Represent a variable list of arguments for piping on a command line, e.g. sam to bam to sorted bam."""

    def __str__(self):
        assert isinstance(self.value, list), \
            "Arguments should be a list"
        assert self.value, "Requires at least one argument"
        # A leading pipe is required so that commands following the last filename
        # do not appear merged.
        # e.g.:  samtools view -bS - | samtools sort -o out.sorted.bam -  [without leading pipe][Incorrect]
        #        | samtools view -bS - | samtools sort -o out.sorted.bam -  [with leading pipe][Correct]
        if any(not isinstance(x, basestring) for x in self.value):
            # Correct for non-string commands.
            # e.g. command classes like Bio.Sequencing.Applications.SamtoolsViewCommandLine
            self.value = map(str, self.value)
        return "| " + " | ".join(self.value)


class _Bowtie2BaseCommandLine(AbstractCommandline):
    """Base bowtie wrapper"""

    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
            _Switch(["-h", "h"],
                    "Print USAGE and DESCRIPTION;  ignore other arguments."),
            _Switch(["--help", "help"],
                    "Print USAGE, DESCRIPTION and ARGUMENTS description; "
                    "ignore other arguments."),
            _Switch(["--version", "version"],
                    "Print version number;  ignore other arguments."),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)

    def _validate(self):
        AbstractCommandline._validate(self)

    def _validate_incompatibilities(self, incompatibles):
        """Used by the bowtie _validate method (PRIVATE)."""
        for element in incompatibles:
            if type(element) is list:
                i = [a for a in element if self._get_parameter(a)]
                if len(i) > 1:
                    raise ValueError("Options {} are incompatible".format(" and ".join(i)))
            elif type(incompatibles) is dict:
                if self._get_parameter(element):
                    for b in incompatibles[element]:
                        if self._get_parameter(b):
                            raise ValueError("Options %s and %s are incompatible."
                                             % (element, b))
            else:
                for a in element:
                    if self._get_parameter(a):
                        for b in incompatibles[a]:
                            if self._get_parameter(b):
                                raise ValueError("Options %s and %s are incompatible."
                                                 % (a, b))


class Bowtie2CommandLine(_Bowtie2BaseCommandLine):
    """Base Bowtie2 wrapper"""

    def __init__(self, cmd='bowtie2', **kwargs):
        assert cmd is not None
        self.parameters = [
            _Option(["-x", "bt2"],
                    "The basename of the index for the reference genome. The basename is the name of any of the index "
                    "files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. bowtie2 looks for the "
                    "specified index first in the current directory, then in the directory specified in the "
                    "BOWTIE2_INDEXES environment variable",
                    filename=True,
                    equate=False),

            _Option(["-1", "m1"],
                    "Comma-separated list of files containing mate 1s (filename usually includes _1), "
                    "e.g. -1 flyA_1.fq,flyB_1.fq. Sequences specified with this option must correspond file-for-file "
                    "and read-for-read with those specified in <m2>. Reads may be a mix of different lengths. If - is "
                    "specified, bowtie2 will read the mate 1s from the standard in or stdin filehandle",
                    equate=False),

            _Option(["-2", "m2"],
                    "Comma-separated list of files containing mate 2s (filename usually includes _2), "
                    "e.g. -2 flyA_2.fq,flyB_2.fq. Sequences specified with this option must correspond file-for-file "
                    "and read-for-read with those specified in <m1>. Reads may be a mix of different lengths. If - is "
                    "specified, bowtie2 will read the mate 2s from the standard in or stdin filehandle",
                    equate=False),

            _Option(["-U", "U"],
                    "Comma-separated list of files containing unpaired reads to be aligned, e.g. lane1.fq,lane2.fq,"
                    "lane3.fq,lane4.fq. Reads may be a mix of different lengths. If - is specified, bowtie2 gets the "
                    "reads from the standard in or stdin filehandle",
                    equate=False),
            _Option(['-S', 'S'],
                    "File to write SAM alignments to. By default, alignments are written to the standard out or "
                    "stdout filehandle (i.e. the console)",
                    filename=True,
                    equate=False)

        ]
        extra_parameters = [
            # Other options
            _Option(["--seed", "seed"],
                    "Use <int> as the seed for pseudo-random number generator. Default: 0",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["--non-deterministic", "non_deterministic"],
                    "Normally, Bowtie 2 re-initializes its pseudo-random generator for each read. It seeds the "
                    "generator with a number derived from (a) the read name, (b) the nucleotide sequence, "
                    "(c) the quality sequence, (d) the value of the --seed option. This means that if two reads are "
                    "identical (same name, same nucleotides, same qualities) Bowtie 2 will find and report the same "
                    "alignment(s) for both, even if there was ambiguity. When --non-deterministic is specified, "
                    "Bowtie 2 re-initializes its pseudo-random generator for each read using the current time. This "
                    "means that Bowtie 2 will not necessarily report the same alignment for two identical reads. This "
                    "is counter-intuitive for some users, but might be more appropriate in situations where the input "
                    "consists of many identical reads"),

            _Switch(["--qc-filter", "qc_filter"],
                    "Filter out reads for which the QSEQ filter field is non-zero. Only has an effect when read "
                    "format is --qseq. Default: off"),

            # Input Options
            _Switch(["-q", "fastq"],
                    "Reads (specified with <m1>, <m2>, <s>) are FASTQ files. FASTQ files usually have "
                    "at. See also: --solexa-quals and --int-quals."),
            _Switch(["--qseq", "qseq"],
                    "Reads (specified with <m1>, <m2>, <s>) are QSEQ files. QSEQ files usually end in s."),
            _Switch(["-f", "fasta"],
                    "Reads (specified with <m1>, <m2>, <s>) are FASTA files. FASTA files usually have "
                    "ore-quals is also set."),
            _Switch(["-r", "unformated"],
                    "Reads (specified with <m1>, <m2>, <s>) are unformated files. With one input sequence per "
                    "if --ignore-quals is also set."),
            _Switch(["-c", "csv"],
                    "The read sequences are given on command line. I.e. <m1>, <m2> and <singles> are CSV files of "
                    "reads rather than lists of or qualities, so -c also implies --ignore-quals."),
            _Switch(["--phred33", "phred33"],
                    "Input qualities are ASCII chars equal to the Phred quality plus 33. This is also called "
                    "the Phred+33 encoding, which is used by the very latest Illumina pipelines"),

            _Switch(["--phred64", "phred64"],
                    "Input qualities are ASCII chars equal to the Phred quality plus 64. This is also called "
                    "the Phred+64 encoding"),

            _Switch(["--solexa-quals", "solexa_quals"],
                    "Convert input qualities from Solexa (which can be negative) to Phred (which can't). This "
                    "scheme was used in older Illumina GA Pipeline versions (prior to 1.3). Default: off"),

            _Switch(["--int-quals", "int_quals"],
                    "Quality values are represented in the read input file as space-separated ASCII integers, "
                    "e.g., 40 40 30 40..., rather than ASCII characters, e.g., II?I.... Integers are treated as "
                    "being on the Phred quality scale unless --solexa-quals is also specified. Default: off"),

            # Preset options in --end-to-end mode
            _Switch(["--very-fast", "very_fast"],
                    "Same as: -D 5 -R 1 -N 0 -L 22 -i S,0,2.50"),

            _Switch(["--fast", "fast"],
                    "Same as: -D 10 -R 2 -N 0 -L 22 -i S,0,2.50"),

            _Switch(["--sensitive", "sensitive"],
                    "Same as: -D 15 -R 2 -L 22 -i S,1,1.15 (default in --end-to-end mode)"),

            _Switch(["--very-sensitive", "very_sensitive"],
                    "Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50"),

            # Preset options in --local mode
            _Switch(["--very-fast-local", "very_fast_local"],
                    "Same as: -D 5 -R 1 -N 0 -L 25 -i S,1,2.00"),

            _Switch(["--fast-local", "fast_local"],
                    "Same as: -D 10 -R 2 -N 0 -L 22 -i S,1,1.75"),

            _Switch(["--sensitive-local", "sensitive_local"],
                    "Same as: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default in --local mode)"),

            _Switch(["--very-sensitive-local", "very_sensitive_local"],
                    "Same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50"),
            # Input configuration options
            _Option(["--skip", "skip"],
                    "Skip (i.e. do not align) the first <int> reads or "
                    "pairs in the input",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--qupto", "qupto"],
                    "Align the first <int> reads or read pairs from the"
                    " input (after the -s/--skip reads or pairs have been skipped), then stop. Default: no limit",
                    checker_function=lambda value: type(value) is int,
                    equate=False),

            _Option(["--trim5", "trim5"],
                    "Trim <int> bases from 5' (left) end of each read before alignment (default: 0)",
                    checker_function=lambda value: type(value) is int,
                    equate=False),

            _Option(["--trim3", "trim3"],
                    "Trim <int> bases from 3' (right) end of each read before alignment (default: 0)",
                    checker_function=lambda value: type(value) is int,
                    equate=False),

            # Alignment options
            _Option(["-N", "num_mismatches"],
                    "Sets the number of mismatches to allowed in a seed alignment during multiseed "
                    "alignment. Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower) "
                    "but increases sensitivity. Default: 0",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["-L", "seed_length"],
                    "Sets the length of the seed substrings to align during multiseed alignment. "
                    "Smaller values make alignment slower but more senstive. Default: the --sensitive preset is used "
                    "by default, which sets -L to 20 both in --end-to-end mode and in --local mode",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["-i", "i_func"],
                    "Sets a function governing the interval between seed substrings to use during multiseed alignment. "
                    "For instance, if the read has 30 characters, and seed length is 10, and the seed interval is 6, "
                    "the seeds extracted will be: Since it's best to use longer intervals for longer reads, this "
                    "parameter sets the interval as a function of the read length, rather than a single one-size-fits-"
                    "all number. For instance, specifying -i S,1,2.5 sets the interval "
                    "function f to f(x) = 1 + 2.5 * sqrt(x), where x is the read length. "
                    "See also: setting function options. If the function returns a result less than 1, it is rounded up"
                    " to 1. Default: the --sensitive preset is used by default, which sets -i to S,1,1.15 "
                    "in --end-to-end mode to -i S,1,0.75 in --local mode.",
                    checker_function=lambda value: re.match('^[CLSG],[-\d\.],[-\d\.]', value) is not None,
                    equate=False),
            _Option(["--n-ceil", "n_ceil"],
                    "Sets a function governing the maximum number of ambiguous characters (usually Ns and/or .s) "
                    "allowed in a read as a function of read length. For instance, specifying -L,0,0.15 sets the "
                    "N-ceiling function f to f(x) = 0 + 0.15 * x, where x is the read length. See also: setting "
                    "function options. Reads exceeding this ceiling are filtered out. Default: L,0,0.15.",
                    checker_function=lambda value: re.match('^[CLSG],[-\d\.],[-\d\.]', value) is not None,
                    equate=False),
            _Option(["--gbar", "gbar"],
                    "Disallow gaps within <int> positions of the beginning or end of the read. Default: 4.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--dpad", "dpad"],
                    "Pads dynamic programming problems by <int> columns on either side to allow gaps. Default: 15.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["--ignore-quals", "ignore_quals"],
                    "When calculating a mismatch penalty, always consider the quality value at the mismatched position "
                    "to be the highest possible, regardless of the actual value. I.e. input is treated as though all "
                    "quality values are high. This is also the default behavior when the input doesn't specify quality "
                    "values (e.g. in -f, -r, or -c modes)"),
            _Switch(["--nofw", "nofw"],
                    "If --nofw is specified, bowtie2 will not attempt to align unpaired reads to the forward (Watson) "
                    "reference strand. In paired-end mode, --nofw and --norc pertain to the fragments; i.e. specifying "
                    "--nofw causes bowtie2 to explore only those paired-end configurations corresponding to fragments "
                    "from the reverse-complement (Crick) strand. Default: both strands enabled"),
            _Switch(["--norc", "norc"],
                    "If --norc is specified, bowtie2 will not attempt to align unpaired reads against the reverse-"
                    "complement Crick reference strand. In paired-end mode, --nofw and --norc pertain to the fragments;"
                    " i.e. specifying --nofw causes bowtie2 to explore only those paired-end configurations "
                    "corresponding to fragments from the reverse-complement (Crick) strand. Default: both strands"),
            _Switch(["--no-1mm-upfront", "no_1mm_upfront"],
                    "By default, Bowtie 2 will attempt to find either an exact or a 1-mismatch end-to-end alignment"
                    " for the read before trying the multiseed heuristic. Such alignments can be found very quickly,"
                    " and many short read alignments have exact or near-exact end-to-end alignments. However, this can "
                    "lead to unexpected alignments when the user also sets options governing the multiseed heuristic, "
                    "like -L and -N. For instance, if the user specifies -N 0 and -L equal to the length of the read, "
                    "the user will be surprised to find 1-mismatch alignments reported. This option prevents Bowtie 2 "
                    "from searching for 1-mismatch end-to-end alignments before using the multiseed heuristic, which "
                    "leads to the expected behavior when combined with options such as -L and -N. This comes at the "
                    "expense of speed"),
            _Switch(["--end-to-end", "end_to_end"],
                    "In this mode, Bowtie 2 requires that the entire read align from one end to the other, without any "
                    "trimming (or soft clipping) of characters from either end. The match bonus --ma always equals 0 in"
                    " this mode, so all alignment scores are less than or equal to 0, and the greatest possible "
                    "alignment score is 0. This is mutually exclusive with --local. --end-to-end is the default mode"),
            _Switch(["--local", "local"],
                    "In this mode, Bowtie 2 does not require that the entire read align from one end to the other. "
                    "Rather, some characters may be omitted (soft clipped) from the ends in order to achieve the "
                    "greatest possible alignment score. The match bonus --ma is used in this mode, and the best "
                    "possible alignment score is equal to the match bonus (--ma) times the length of the read. "
                    "Specifying --local and one of the presets (e.g. --local --very-fast) is equivalent to specifying "
                    "the local version of the preset (--very-fast-local). This is mutually exclusive with --end-to-end."
                    " --end-to-end is the default mode"),

            # Scoring Options
            _Option(["--score-min", "score_min"],
                    "Sets a function governing the minimum alignment score needed for an alignment to be considered "
                    "valid (i.e. good enough to report). This is a function of read length. For instance, specifying "
                    "L,0,-0.6 sets the minimum-score function f to f(x) = 0 + -0.6 * x, where x is the read length."
                    " See also: setting function options. The default in --end-to-end mode is L,-0.6,-0.6 "
                    "and the default in --local mode is G,20,8.",
                    checker_function=lambda value: re.match('^[CLSG],[-\d\.],[-\d\.]', value) is not None,
                    equate=False),
            _Option(["--ma", "ma"],
                    "Sets the match bonus. In --local mode <int> is added to the alignment score for each "
                    "position where a read character aligns to a reference character and the characters match. "
                    "Not used in --end-to-end mode. Default: 2.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--np", "np"],
                    "Sets penalty for positions where the read, reference, or both, contain an ambiguous "
                    "character such as N. Default: 1.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--rdg", "rdg"],
                    "Sets the read gap open (<int1>) and extend (<int2>) penalties. A read gap of length N gets"
                    " a penalty of <int1> + N * <int2>. Default: 5, 3.",
                    checker_function=lambda value: re.match('[-d.],[-d.]', value) is not None,
                    equate=False),
            _Option(["--rfg", "rfg"],
                    "Sets the reference gap open (<int1>) and extend (<int2>) penalties. A reference gap of "
                    "length N gets a penalty of <int1> + N * <int2>. Default: 5, 3.",
                    checker_function=lambda value: re.match('[-d.],[-d.]', value) is not None,
                    equate=False),
            _Option(["--mp", "mp"],
                    "Sets the maximum (MX) and minimum (MN) mismatch penalties, both integers. A number less "
                    "than or equal to MX and greater than or equal to MN is subtracted from the alignment score for "
                    "each position where a read character aligns to a reference character, the characters do not match,"
                    " and neither is an N. If --ignore-quals is specified, the number subtracted quals MX. "
                    "Otherwise, the number subtracted is MN + floor( (MX-MN)(MIN(Q, 40.0)/40.0) ) "
                    "where Q is the Phred quality value. Default: MX = 6, MN = 2.",
                    checker_function=lambda value: re.match('[-d.],[-d.]', value) is not None,
                    equate=False),

            # Reporting Options
            _Option(["-k", "k"],
                    "By default, bowtie2 searches for distinct, valid alignments for each read. When it finds a"
                    " valid alignment, it continues looking for alignments that are nearly as good or better. The best "
                    "alignment found is reported (randomly selected from among best if tied). Information about the "
                    "best alignments is used to estimate mapping quality and to set SAM optional fields, such as "
                    "AS:i and XS:i.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["-a", "a"],
                    "Like -k but with no upper limit on number of alignments to search for. "
                    "-a is mutually exclusive with -k."),

            # Effort Options
            _Option(["-D", "D"],
                    "Up to <int> consecutive seed extension attempts can fail before Bowtie 2 moves on, using"
                    " the alignments found so far. A seed extension fails if it does not yield a new best or a new "
                    "second-best alignment. This limit is automatically adjusted up when -k or -a are specified. "
                    "Default: 15.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["-R", "R"],
                    "<int> is the maximum number of times Bowtie 2 will re-seed reads with repetitive seeds. "
                    "When re-seeding, Bowtie 2 simply chooses a new set of reads (same length, same number of "
                    "mismatches allowed) at different offsets and searches for more alignments. A read is considered "
                    "to have repetitive seeds if the total number of seed hits divided by the number of seeds that "
                    "aligned at least once is greater than 300. Default: 2.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),

            # Paired-end options
            _Option(["--minins", "minins"],
                    "The minimum fragment length for valid paired-end alignments. E.g. if -I 60 is specified "
                    "and a paired-end alignment consists of two 20-bp alignments in the appropriate orientation with "
                    "a 20-bp gap between them, that alignment is considered valid (as long as -X is also satisfied). "
                    "A 19-bp gap would not be valid in that case. If trimming options -3 or -5 are also used, "
                    "the -I constraint is applied with respect to the untrimmed mates. The larger the difference "
                    "between -I and -X, the slower Bowtie 2 will run. This is because larger differences bewteen -I "
                    "and -X require that Bowtie 2 scan a larger window to determine if a concordant alignment exists. "
                    "For typical fragment length ranges (200 to 400 nucleotides), Bowtie 2 is very efficient. "
                    "Default: 0 (essentially imposing no minimum)",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--maxins", "maxins"],
                    "The maximum fragment length for valid paired-end alignments. E.g. if -X 100 is specified "
                    "and a paired-end alignment consists of two 20-bp alignments in the proper orientation with a "
                    "60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied). "
                    "A 61-bp gap would not be valid in that case. If trimming options -3 or -5 are also used, the "
                    "-X constraint is applied with respect to the untrimmed mates, not the trimmed mates. The larger "
                    "the difference between -I and -X, the slower Bowtie 2 will run. This is because larger differences"
                    " bewteen -I and -X require that Bowtie 2 scan a larger window to determine if a concordant "
                    "alignment exists. For typical fragment length ranges (200 to 400 nucleotides), "
                    "Bowtie 2 is very efficient. Default: 500",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["--fr", "fr"],
                    "The upstream/downstream mate orientations for a valid paired-end alignment against the "
                    "forward reference strand. E.g., if --fr is specified and there is a candidate paired-end "
                    "alignment where mate 1 appears upstream of the reverse complement of mate 2 and the fragment "
                    "length constraints (-I and -X) are met, that alignment is valid. Also, if mate 2 appears "
                    "upstream of the reverse complement of mate 1 and all other constraints are met, "
                    "that too is valid. --rf likewise requires that an upstream mate1 be reverse-complemented "
                    "and a downstream mate2 be forward-oriented. --ff requires both an upstream mate 1 and a "
                    "downstream mate 2 to be forward-oriented. "
                    "Default: --fr (appropriate for Illumina's Paired-end Sequencing Assay)."),

            _Switch(["--rf", "rf"],
                    "The upstream/downstream mate orientations for a valid paired-end alignment against the "
                    "forward reference strand. E.g., if --fr is specified and there is a candidate paired-end "
                    "alignment where mate 1 appears upstream of the reverse complement of mate 2 and the fragment "
                    "length constraints (-I and -X) are met, that alignment is valid. Also, if mate 2 appears "
                    "upstream of the reverse complement of mate 1 and all other constraints are met, "
                    "that too is valid. --rf likewise requires that an upstream mate1 be reverse-complemented "
                    "and a downstream mate2 be forward-oriented. --ff requires both an upstream mate 1 and a "
                    "downstream mate 2 to be forward-oriented. "
                    "Default: --fr (appropriate for Illumina's Paired-end Sequencing Assay)."),
            _Switch(["--ff", "ff"],
                    "The upstream/downstream mate orientations for a valid paired-end alignment against the "
                    "forward reference strand. E.g., if --fr is specified and there is a candidate paired-end "
                    "alignment where mate 1 appears upstream of the reverse complement of mate 2 and the fragment "
                    "length constraints (-I and -X) are met, that alignment is valid. Also, if mate 2 appears "
                    "upstream of the reverse complement of mate 1 and all other constraints are met, "
                    "that too is valid. --rf likewise requires that an upstream mate1 be reverse-complemented "
                    "and a downstream mate2 be forward-oriented. --ff requires both an upstream mate 1 and a "
                    "downstream mate 2 to be forward-oriented. "
                    "Default: --fr (appropriate for Illumina's Paired-end Sequencing Assay)."),
            _Switch(["--no-mixed", "no_mixed"],
                    "By default, when bowtie2 cannot find a concordant or discordant alignment for a pair, it "
                    "then tries to find alignments for the individual mates. This option disables that behavior."),

            _Switch(["--no-discordant", "no_discordant"],
                    "By default, bowtie2 looks for discordant alignments if it cannot find any concordant "
                    "alignments. A discordant alignment is an alignment where both mates align uniquely, "
                    "but that does not satisfy the paired-end constraints (--fr/--rf/--ff, -I, -X). "
                    "This option disables that behavior."),

            _Switch(["--dovetail", "dovetail"],
                    "If the mates dovetail, that is if one mate alignment extends past the beginning of the "
                    "other such that the wrong mate begins upstream, consider that to be concordant. See also: "
                    "Mates can overlap, contain or dovetail each other. Default: mates cannot dovetail "
                    "in a concordant alignment."),

            _Switch(["--no-contain", "no_contain"],
                    "If one mate alignment contains the other, consider that to be non-concordant. See also: "
                    "Mates can overlap, contain or dovetail each other. Default: a mate can contain "
                    "the other in a concordant alignment."),

            _Switch(["--no-overlap", "no_overlap"],
                    "If one mate alignment overlaps the other at all, consider that to be non-concordant. See "
                    "also: Mates can overlap, contain or dovetail each other. Default: mates can overlap in "
                    "a concordant alignment."),

            # SAM options
            _Switch(["--no-unal", "no_unal"],
                    "Suppress SAM records for reads that failed to align"),
            _Switch(["--no-hd", "no_hd"],
                    "Suppress SAM header lines (starting with"),
            _Switch(["--no-sq", "no_sq"],
                    "Suppress @SQ SAM header lines"),
            _Switch(["--omit-sec-seq", "omit_sec_seq"],
                    "When printing secondary alignments, Bowtie 2 by default will write out the SEQ and QUAL strings. "
                    "Specifying this option causes Bowtie 2 to print an asterix in those fields instead."),
            _Option(["--rg-id", "rg_id"],
                    "Set the read group ID to <text>. This causes the SAM @RG header line to be printed, with <text> as"
                    " the value associated with the ID: tag. It also causes the RG:Z: extra field to be attached to "
                    "each SAM output record, with value set to <text>.",
                    checker_function=lambda value: type(value) is str,
                    equate=False),
            _Option(["--rg", "rg"],
                    "Add <text> (usually of the form TAG:VAL, e.g. SM:Pool1) as a field on the @RG header line. Note: "
                    "in order for the @RG line to appear, --rg-id must also be specified. This is because the ID tag is"
                    " required by the SAM Spec. Specify --rg multiple times to set multiple fields. See the SAM "
                    "Spec for details about what fields are legal.",
                    checker_function=lambda value: type(value) is str,
                    equate=False),

            # Output options
            _Option(["--un", "un"],
                    "Write unpaired reads that fail to align to file at <path>. These reads correspond to the SAM "
                    "records with the FLAGS 0x4 bit set and neither the 0x40 nor 0x80 bits set. Reads written in this "
                    "way will appear exactly as they did in the input file, without any modification (same sequence, "
                    "same name, same quality string, same quality encoding). Reads will not necessarily appear in the "
                    "same order as they did in the input",
                    filename=True,
                    equate=False),
            _Option(["--un-gz", "un_gz"],
                    "Write unpaired reads that fail to align to file at <path>. These reads correspond to the SAM "
                    "records with the FLAGS 0x4 bit set and neither the 0x40 nor 0x80 bits set. If --un-gz is "
                    "specified, output will be gzip compressed. Reads written in this way will appear exactly as they "
                    "did in the input file, without any modification (same sequence, same name, same quality string, "
                    "same quality encoding). Reads will not necessarily appear in the same order as they did in the "
                    "input",
                    filename=True,
                    equate=False),
            _Option(["--un-bz2", "un_bz2"],
                    "Write unpaired reads that fail to align to file at <path>. These reads correspond to the SAM "
                    "records with the FLAGS 0x4 bit set and neither the 0x40 nor 0x80 bits set.  If --un-bz2 is "
                    "specified, output will be bzip2 compressed. Reads written in this way will appear exactly as "
                    "they did in the input file, without any modification (same sequence, same name, same quality "
                    "string, same quality encoding). Reads will not necessarily appear in the same order as they did "
                    "in the input",
                    filename=True,
                    equate=False),
            _Option(["--un-lz4", "un_lz4"],
                    "Write unpaired reads that fail to align to file at <path>. These reads correspond to the SAM "
                    "records with the FLAGS 0x4 bit set and neither the 0x40 nor 0x80 bits set. If --un-lz4 is "
                    "specified, output will be lz4 compressed. Reads written in this way will appear exactly as they "
                    "did in the input file, without any modification (same sequence, same name, same quality string, "
                    "same quality encoding). Reads will not necessarily appear in the same order as they did in the "
                    "input",
                    filename=True,
                    equate=False),

            _Option(["--al", "al"],
                    "Write unpaired reads that align at least once to file at <path>. These reads correspond to the "
                    "SAM records with the FLAGS 0x4, 0x40, and 0x80 bits unset. Reads written in this way will appear "
                    "exactly as they did in the input file, without any modification (same sequence, same name, "
                    "same quality string, same quality encoding). Reads will not necessarily appear in the same order "
                    "as they did in the input",
                    filename=True,
                    equate=False),
            _Option(["--al-gz", "al_gz"],
                    "Write unpaired reads that align at least once to file at <path>. These reads correspond to the "
                    "SAM records with the FLAGS 0x4, 0x40, and 0x80 bits unset. If --al-gz is specified, output will "
                    "be gzip compressed. Reads written in this way will appear exactly as they did in the input file, "
                    "without any modification (same sequence, same name, same quality string, same quality encoding). "
                    "Reads will not necessarily appear in the same order as they did in the input",
                    filename=True,
                    equate=False),
            _Option(["--al-bz2", "al_bz2"],
                    "Write unpaired reads that align at least once to file at <path>. These reads correspond to the "
                    "SAM records with the FLAGS 0x4, 0x40, and 0x80 bits unset. If --al-bz2 is specified, output will "
                    "be bzip2 compressed. Reads written in this way will appear exactly as they did in the input "
                    "file, without any modification (same sequence, same name, same quality string, same quality "
                    "encoding). Reads will not necessarily appear in the same order as they did in the input",
                    filename=True,
                    equate=False),
            _Option(["--al-lz4", "al_lz4"],
                    "Write unpaired reads that align at least once to file at <path>. These reads correspond to the "
                    "SAM records with the FLAGS 0x4, 0x40, and 0x80 bits unset. If --al-lz4 is specified, output will "
                    "be lz4 compressed. Reads written in this way will appear exactly as they did in the input file, "
                    "without any modification (same sequence, same name, same quality string, same quality encoding). "
                    "Reads will not necessarily appear in the same order as they did in the input",
                    filename=True,
                    equate=False),

            _Option(["--un-conc", "un_conc"],
                    "Write paired-end reads that fail to align concordantly to file(s) at <path>. These reads "
                    "correspond to the SAM records with the FLAGS 0x4 bit set and either the 0x40 or 0x80 bit set ("
                    "depending on whether it's mate #1 or #2). .1 and .2 strings are added to the filename to "
                    "distinguish which file contains mate #1 and mate #2. If a percent symbol, %, is used in <path>, "
                    "the percent symbol is replaced with 1 or 2 to make the per-mate filenames. Otherwise, "
                    ".1 or .2 are added before the final dot in <path> to make the per-mate filenames. Reads written "
                    "in this way will appear exactly as they did in the input files, without any modification (same "
                    "sequence, same name, same quality string, same quality encoding). Reads will not necessarily "
                    "appear in the same order as they did in the inputs",
                    filename=True,
                    equate=False),
            _Option(["--un-conc-gz", "un_conc_gz"],
                    "Write paired-end reads that fail to align concordantly to file(s) at <path>. These reads "
                    "correspond to the SAM records with the FLAGS 0x4 bit set and either the 0x40 or 0x80 bit set ("
                    "depending on whether it's mate #1 or #2). .1 and .2 strings are added to the filename to "
                    "distinguish which file contains mate #1 and mate #2. If a percent symbol, %, is used in <path>, "
                    "the percent symbol is replaced with 1 or 2 to make the per-mate filenames. Otherwise, "
                    ".1 or .2 are added before the final dot in <path> to make the per-mate filenames. Reads written "
                    "in this way will appear exactly as they did in the input files, without any modification (same "
                    "sequence, same name, same quality string, same quality encoding). Reads will not necessarily "
                    "appear in the same order as they did in the inputs",
                    filename=True,
                    equate=False),
            _Option(["--un-conc-bz2", "un_conc_bz2"],
                    "Write paired-end reads that fail to align concordantly to file(s) at <path>. These reads "
                    "correspond to the SAM records with the FLAGS 0x4 bit set and either the 0x40 or 0x80 bit set ("
                    "depending on whether it's mate #1 or #2). .1 and .2 strings are added to the filename to "
                    "distinguish which file contains mate #1 and mate #2. If a percent symbol, %, is used in <path>, "
                    "the percent symbol is replaced with 1 or 2 to make the per-mate filenames. Otherwise, "
                    ".1 or .2 are added before the final dot in <path> to make the per-mate filenames. Reads written "
                    "in this way will appear exactly as they did in the input files, without any modification (same "
                    "sequence, same name, same quality string, same quality encoding). Reads will not necessarily "
                    "appear in the same order as they did in the inputs",
                    filename=True,
                    equate=False),
            _Option(["--un-conc-lz4", "un_conc_lz4"],
                    "Write paired-end reads that fail to align concordantly to file(s) at <path>. These reads "
                    "correspond to the SAM records with the FLAGS 0x4 bit set and either the 0x40 or 0x80 bit set ("
                    "depending on whether it's mate #1 or #2). .1 and .2 strings are added to the filename to "
                    "distinguish which file contains mate #1 and mate #2. If a percent symbol, %, is used in <path>, "
                    "the percent symbol is replaced with 1 or 2 to make the per-mate filenames. Otherwise, "
                    ".1 or .2 are added before the final dot in <path> to make the per-mate filenames. Reads written "
                    "in this way will appear exactly as they did in the input files, without any modification (same "
                    "sequence, same name, same quality string, same quality encoding). Reads will not necessarily "
                    "appear in the same order as they did in the inputs",
                    filename=True,
                    equate=False),

            _Option(["--al-conc", "al_conc"],
                    "Write paired-end reads that align concordantly at least once to file(s) at <path>. These reads "
                    "correspond to the SAM records with the FLAGS 0x4 bit unset and either the 0x40 or 0x80 bit set ("
                    "depending on whether it's mate #1 or #2). .1 and .2 strings are added to the filename to "
                    "distinguish which file contains mate #1 and mate #2. If a percent symbol, %, is used in <path>, "
                    "the percent symbol is replaced with 1 or 2 to make the per-mate filenames. Otherwise, "
                    ".1 or .2 are added before the final dot in <path> to make the per-mate filenames. Reads written "
                    "in this way will appear exactly as they did in the input files, without any modification (same "
                    "sequence, same name, same quality string, same quality encoding). Reads will not necessarily "
                    "appear in the same order as they did in the inputs",
                    filename=True,
                    equate=False),
            _Option(["--al-conc-gz", "al_conc_gz"],
                    "Write paired-end reads that align concordantly at least once to file(s) at <path>. These reads "
                    "correspond to the SAM records with the FLAGS 0x4 bit unset and either the 0x40 or 0x80 bit set ("
                    "depending on whether it's mate #1 or #2). .1 and .2 strings are added to the filename to "
                    "distinguish which file contains mate #1 and mate #2. If a percent symbol, %, is used in <path>, "
                    "the percent symbol is replaced with 1 or 2 to make the per-mate filenames. Otherwise, "
                    ".1 or .2 are added before the final dot in <path> to make the per-mate filenames. Reads written "
                    "in this way will appear exactly as they did in the input files, without any modification (same "
                    "sequence, same name, same quality string, same quality encoding). Reads will not necessarily "
                    "appear in the same order as they did in the inputs",
                    filename=True,
                    equate=False),
            _Option(["--al-conc-bz2", "al_conc_bz2"],
                    "Write paired-end reads that align concordantly at least once to file(s) at <path>. These reads "
                    "correspond to the SAM records with the FLAGS 0x4 bit unset and either the 0x40 or 0x80 bit set ("
                    "depending on whether it's mate #1 or #2). .1 and .2 strings are added to the filename to "
                    "distinguish which file contains mate #1 and mate #2. If a percent symbol, %, is used in <path>, "
                    "the percent symbol is replaced with 1 or 2 to make the per-mate filenames. Otherwise, "
                    ".1 or .2 are added before the final dot in <path> to make the per-mate filenames. Reads written "
                    "in this way will appear exactly as they did in the input files, without any modification (same "
                    "sequence, same name, same quality string, same quality encoding). Reads will not necessarily "
                    "appear in the same order as they did in the inputs",
                    filename=True,
                    equate=False),
            _Option(["--al-conc-lz4", "al_conc_lz4"],
                    "Write paired-end reads that align concordantly at least once to file(s) at <path>. These reads "
                    "correspond to the SAM records with the FLAGS 0x4 bit unset and either the 0x40 or 0x80 bit set ("
                    "depending on whether it's mate #1 or #2). .1 and .2 strings are added to the filename to "
                    "distinguish which file contains mate #1 and mate #2. If a percent symbol, %, is used in <path>, "
                    "the percent symbol is replaced with 1 or 2 to make the per-mate filenames. Otherwise, "
                    ".1 or .2 are added before the final dot in <path> to make the per-mate filenames. Reads written "
                    "in this way will appear exactly as they did in the input files, without any modification (same "
                    "sequence, same name, same quality string, same quality encoding). Reads will not necessarily "
                    "appear in the same order as they did in the inputs",
                    filename=True,
                    equate=False),

            _Option(["--met-file", "met_file"],
                    "Write bowtie2 metrics to file <path>. Having alignment metric can be useful for debugging "
                    "certain problems, especially performance issues. See also: --met. Default: metrics disabled",
                    filename=True,
                    equate=False),

            _Option(["--met-stderr", "met_stderr"],
                    "Write bowtie2 metrics to the standard error (stderr) filehandle. This is not mutually exclusive "
                    "with --met-file. Having alignment metric can be useful for debugging certain problems, "
                    "especially performance issues. See also: --met. Default: metrics disabled",
                    filename=True,
                    equate=False),

            _Option(["--met", "met"],
                    "Write a new bowtie2 metrics record every <int> seconds. Only matters if either --met-stderr or "
                    "--met-file are specified. Default: 1.",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["--time", "time"],
                    "Print the wall-clock time required to load the index files and align the reads. This is printed "
                    "to the standard error (stderr) filehandle. Default: off"),
            _Switch(["--quiet", "quiet"],
                    "Print nothing besides alignments and serious errors"),
            # Preformance Options
            _Option(["--offrate", "offrate"],
                    "Override the offrate of the index with <int>. If <int> is greater than the offrate used to build "
                    "the index, then some row markings are discarded when the index is read into memory. This reduces "
                    "the memory footprint of the aligner but requires more time to calculate text offsets. <int> must "
                    "be greater than the value used to build the index",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--threads", "threads"],
                    "Launch NTHREADS parallel search threads (default: 1). Threads will run on separate "
                    "processors/cores and synchronize when parsing reads and outputting alignments. Searching for "
                    "alignments is highly parallel, and speedup is close to linear. Increasing -p increases Bowtie "
                    "2's memory footprint. E.g. when aligning to a human genome index, increasing -p from 1 to 8 "
                    "increases the memory footprint by a few hundred megabytes. This option is only available if "
                    "bowtie is linked with the pthreads library (i.e. if BOWTIE_PTHREADS=0 is not specified at build "
                    "time)",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["--reorder", "reorder"],
                    "Guarantees that output SAM records are printed in an order corresponding to the order of the "
                    "reads in the original input file, even when -p is set greater than 1. Specifying --reorder and "
                    "setting -p greater than 1 causes Bowtie 2 to run somewhat slower and use somewhat more memory "
                    "then if --reorder were not specified. Has no effect if -p is set to 1, since output order will "
                    "naturally correspond to input order in that case"),
            _Switch(["--mm", "mm"],
                    "Use memory-mapped I/O to load the index, rather than typical file I/O. Memory-mapping allows "
                    "many concurrent bowtie processes on the same computer to share the same memory image of the "
                    "index (i.e. you pay the memory overhead just once). This facilitates memory-efficient "
                    "parallelization of bowtie in situations where using -p is not possible or not preferable"),
        ]
        pipe_parameters = [
            _PipeArgumentList(["samcmds", "samtools"],
                              "Allow user to pipe bowtie2 output to samtools for bam output")

        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            # add pipe parameters to the end
            self.parameters = extra_parameters + self.parameters + pipe_parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters + pipe_parameters
        _Bowtie2BaseCommandLine.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = [["local", "end_to_end"],
                         ["k", "a"],
                         ["al", "al_gz", "al_bz2", "al_lz4"],
                         ["un", "un_gz", "un_bz2", "un_lz4"],
                         ["un_conc", "un_conc_gz", "un_conc_bz2", "un_conc_lz4"],
                         ["al_conc", "al_conc_gz", "al_conc_bz2", "al_lz4"]]
        self._validate_incompatibilities(incompatibles)
        # TODO add incompatibilites
        if self.bt2:
            if (not self.m1 and not self.m2) and not self.U:
                raise ValueError("Option bowtie2 requires input fastq.")
        _Bowtie2BaseCommandLine._validate(self)


class _Bowtie2SeqBaseCommandLine(_Bowtie2BaseCommandLine):
    """Base bowtie wrapper"""

    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        self.parameters += [
            _Argument(["bt2"],
                      "bt2 filename minus trailing .1.bt2/.2.bt2. bt2 data to files with this dir/basename")
        ]
        extra_parameters = [
            _Switch(["--large-index", "large_index"],
                    "Force bowtie2-build to build a large index, even if the reference is less than ~ 4 billion "
                    "nucleotides inlong."),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _Bowtie2BaseCommandLine.__init__(self, cmd, **kwargs)

    def _validate(self):
        incompatibles = []
        self._validate_incompatibilities(incompatibles)
        _Bowtie2BaseCommandLine._validate(self)


class Bowtie2BuildCommandLine(_Bowtie2SeqBaseCommandLine):
    """Base bowtie2-build wrapper"""

    def __init__(self, cmd='bowtie2-build', **kwargs):
        assert cmd is not None
        self.parameters = [
            _Argument(["reference"],
                      "comma-separated list of files with ref sequences")
        ]
        extra_parameters = [
            # Other options
            _Option(["--seed", "seed"],
                    "Use <int> as the seed for pseudo-random number generator. Default: 0",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["--quiet", "quiet"],
                    "Print nothing besides alignments and serious errors"),
            _Option(["--bmax", "bmax"],
                    "The maximum number of suffixes allowed in a block. Allowing more suffixes per block makes "
                    "indexing faster, but increases peak memory usage. Setting this option overrides any previous "
                    "setting for --bmax, or --bmaxdivn. Default (in terms of the --bmaxdivn parameter) is --bmaxdivn "
                    "4. This is configured automatically by default; use -a/--noauto to configure manually",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--bmaxdivn", "bmaxdivn"],
                    "The maximum number of suffixes allowed in a block, expressed as a fraction of the length of the "
                    "reference. Setting this option overrides any previous setting for --bmax, or --bmaxdivn. "
                    "Default: --bmaxdivn 4. This is configured automatically by default; use -a/--noauto to configure "
                    "manually",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--dcv", "dcv"],
                    "Use <int> as the period for the difference-cover sample. A larger period yields less memory "
                    "overhead, but may make suffix sorting slower, especially if repeats are present. Must be a power "
                    "of 2 no greater than 4096. Default: 1024. This is configured automatically by default; use "
                    "-a/--noauto to configure manually",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--offrate", "offrate"],
                    "To map alignments back to positions on the reference sequences, it's necessary to annotate ("
                    "mark) some or all of the Burrows-Wheeler rows with their corresponding location on the genome. "
                    "-o/--offrate governs how many rows get marked: the indexer will mark every 2^<int> rows. Marking "
                    "more rows makes reference-position lookups faster, but requires more memory to hold the "
                    "annotations at runtime. The default is 5 (every 32nd row is marked; for human genome, "
                    "annotations occupy about 340 megabytes)",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--ftabchars", "ftabchars"],
                    "The ftab is the lookup table used to calculate an initial Burrows-Wheeler range with respect to "
                    "the first <int> characters of the query. A larger <int> yields a larger lookup table but faster "
                    "query times. The ftab has size 4^(<int>+1) bytes. The default setting is 10 (ftab is 4MB)",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Option(["--cutoff", "cutoff"],
                    "Index only the first <int> bases of the reference sequences (cumulative across sequences) and "
                    "ignore the rest",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["-f", "f"],
                    "The reference input files (specified as <reference_in>) are FASTA files (usually having "
                    "extension .fa, .mfa, .fna or similar)."),
            _Switch(["-c", "c"],
                    "The reference sequences are given on the command line. I.e. <reference_in> is a comma-separated "
                    "list of sequences rather than a list of FASTA files."),
            _Switch(["--noauto", "noauto"],
                    "Disable the default behavior whereby bowtie2-build automatically selects values for the --bmax, "
                    "--dcv and --packed parameters according to available memory. Instead, user may specify values "
                    "for those parameters. If memory is exhausted during indexing, an error message will be printed; "
                    "it is up to the user to try new parameters."),
            _Switch(["--packed", "packed"],
                    "Use a packed (2-bits-per-nucleotide) representation for DNA strings. This saves memory but makes "
                    "indexing 2-3 times slower. Default: off. This is configured automatically by default; use "
                    "-a/--noauto to configure manually."),
            _Switch(["--nodc", "nodc"],
                    "Disable use of the difference-cover sample. Suffix sorting becomes quadratic-time in the worst "
                    "case (where the worst case is an extremely repetitive reference). Default: off."),
            _Switch(["--noref", "noref"],
                    "Do not build the NAME.3.bt2 and NAME.4.bt2 portions of the index, which contain a bitpacked "
                    "version of the reference sequences and are used for paired-end alignment."),
            _Switch(["--justref", "justref"],
                    "Build only the NAME.3.bt2 and NAME.4.bt2 portions of the index, which contain a bitpacked "
                    "version of the reference sequences and are used for paired-end alignment."),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _Bowtie2SeqBaseCommandLine.__init__(self, cmd, **kwargs)


class Bowtie2InspectCommandLine(_Bowtie2SeqBaseCommandLine):
    """Base bowtie2-inspoect wrapper"""

    def __init__(self, cmd='bowtie2-inspect', **kwargs):
        assert cmd is not None
        self.parameters = list()
        extra_parameters = [
            _Option(["--across", "across"],
                    "When printing FASTA output, output a newline character every <int> bases (default: 60).",
                    checker_function=lambda value: type(value) is int,
                    equate=False),
            _Switch(["--names", "names"],
                    "Print reference sequence names, one per line, and quit."),
            _Switch(["--summary", "summary"],
                    "Print a summary that includes information about index settings, as well as the names and lengths "
                    "of the input sequences. Fields are separated by tabs. Colorspace is always set to 0 for Bowtie "
                    "2."),
            _Switch(["--verbose", "verbose"],
                    "Print verbose output (for debugging)."),
        ]
        try:
            # Insert extra parameters - at the start just in case there
            # are any arguments which must come last:
            self.parameters = extra_parameters + self.parameters
        except AttributeError:
            # Should we raise an error?  The subclass should have set this up!
            self.parameters = extra_parameters
        _Bowtie2SeqBaseCommandLine.__init__(self, cmd, **kwargs)


if __name__ == '__main__':
    from Bio.Sequencing.Applications import SamtoolsViewCommandline, SamtoolsSortCommandline

    ubam = "/data/2015-SEQ-1283/qualimap_results/2015-SEQ-1283.sorted.bam"
    samsortt = SamtoolsSortCommandline(input_bam="-", out_prefix=ubam[:-4])
    samtoolss = [SamtoolsViewCommandline(b=True, S=True, input_file="-"), samsortt]
    # print samtools
    print Bowtie2CommandLine(bt2="test", m1="none", m2="yes", samtools=samtoolss)
    # print Bowtie2InspectCommandLine(bt2="test")
    pass