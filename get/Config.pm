=head1 LICENSE

Copyright (c) 2007-2011 Illumina, Inc.

This software is covered by the "Illumina Genome Analyzer Software 
License Agreement" and the "Illumina Source Code License Agreement", 
and certain third party copyright/licenses, and any user of this 
source file is bound by the terms therein (see accompanying files
Illumina_Genome_Analyzer_Software_License_Agreement.pdf and
Illumina_Source_Code_License_Agreement.pdf and third party 
copyright/license notices).

This file is part of the Consensus Assessment of Sequence And VAriation
(CASAVA) software package.

=head1 NAME

Casava::Config - Utility library for the management of the config files

=head1 SYNOPSIS

  use Casava::Alignment::Config;
  my $config = Casava::Alignment::Config->new;
  $config->readTxt("/path/to/input/config.txt");
  $config->writeTxt("/path/to/output/config.txt");
  $config->writeXml("/path/to/output/config.xml");
  $config->writeMakefile("/path/to/output/Makefile.config");

=head1 DESCRIPTION

This script is a rewrite of the original GERALD library originally
implemented for the Pipeline.

Library for the management of the configuration files for the
Alignment folder.

=head1 TEXT FILE FORMAT

The text file format include one variable assignment per line. The variable assignment can be global, ro a specific barcode, sample, reference or project. A variable assignment can be prefixed by a lane specification, in which case the assignment applies only to the specified lanes. Lines starting with a '#' are comments.

  # A global assignment
  <VARIABLE_NAME> <variable-value>
  # An assignment for project <PROJECT-NAME>
  PROJECT <PROJECT-NAME> <VARIABLE_NAME> <variable-name>
  # An assignment for sample <SAMPLE_NAME>
  SAMPLE <SAMPLE_NAME> <VARIABLE_NAME> <variable-value>
  # An assignment for reference <REFERENCE_NAME>
  REFERENCE <REFERENCE_NAME> <VARIABLE_NAME> <variable-value>
  # An assignment for barcode <BARCODE_SEQUENCE>
  BARCODE <BARCODE_SEQUENCE> <VARIABLE_NAME> <variable-value>
  # An assignment for lanes 1, 3 and 6 only
  136:<VARIABLE_NAME> <variable-value>
  # An assignment for project <PROJECT-NAME> for lanes 1, 3 and 6 only
  136:PROJECT <PROJECT-NAME> <VARIABLE_NAME> <variable-name>
  # An assignment for sample <SAMPLE_NAME> for lanes 1, 3 and 6 only
  136:SAMPLE <SAMPLE_NAME> <VARIABLE_NAME> <variable-value>
  # An assignment for reference <REFERENCE_NAME> for lanes 1, 3 and 6 only
  136:REFERENCE <REFERENCE_NAME> <VARIABLE_NAME> <variable-value>
  # An assignment for barcode <BARCODE_SEQUENCE> for lanes 1, 3 and 6 only
  136:BARCODE <BARCODE_SEQUENCE> <VARIABLE_NAME> <variable-value>

=head1 XML FILE FORMAT

The xml file format expands the text file format. It has the following entries:

  <RunParameters>
    <FlowCell>
      <Default>
      </Default>
      <Projects>
        <Project name="projectA">
        </Project>
      </Projects>
      <References>
        <Reference name="human36">
        </Reference>
      </References>
      <Samples>
        <Sample name="someSample">
        </Sample>
      </Samples>
      <Barcodes>
        <Barcode name="ACGTAC">
        </Barcode>
      </Barcodes>
    </FlowCell>
    <Lanes>
      <Lane index="1">
        <Default>
        </Default>
        <Projects>
          <Project name="projectA">
          </Project>
        </Projects>
        <References>
          <Reference name="human36">
          </Reference>
        </References>
        <Samples>
          <Sample name="someSample">
          </Sample>
        </Samples>
        <Barcodes>
          <Barcode name="ACGTAC">
          </Barcode>
        </Barcodes>        
      </Lane>
    </Lanes>
  </RunParameters>

=head1 SUBROUTINES

=cut

# POD for for each subroutines is right before the implementation
# More POD after __END__

package Casava::Alignment::Config;

use strict;
use warnings "all";
use Carp;
use XML::Simple;
use File::Spec;
use List::Util qw(min);

use Casava::Common::Log;
use Casava::Common::Utils qw(expandUseBases getFlowCellId);

use Exporter qw( import );
our @EXPORT_OK = qw( new readTxt writeXml toXml getDefaultVariables setDefaultVariable validate validateDataset);

=pod

=head2 new

Create a new instance of the alignment configuration.

I<Parameters:>

=over 4

None

=back

I<Returns:>

An instance to the config component.

I<Exceptions:>

=over 4

None

=back

=cut

sub new
{
    my $class = shift;
    my $self = {Defaults => {},
                Projects => {Project => []},
                References => {Reference => []},
                Samples => {Sample => []},
                Barcodes => {Barcode => []},
                Lanes => {Lane => []},
    };
    bless $self, $class;
    return $self;
}

sub software($;$)
{
    my $self = shift;
    if (@_) { $self->{Software} = shift }
    return $self->{Software};
}

sub cmdAndArgs($;$)
{
    my $self = shift;
    if (@_) { $self->{Software}->{CmdAndArgs} = shift }
    return $self->{Software}->{CmdAndArgs};
}

sub demuxSoftware($;$)
{
    my $self = shift;
    if (@_) { $self->{Software}->{Software} = shift }
    return $self->{Software}->{Software};
}

=head2 readTxt

Load configuration information from a text file.

I<Parameters:>

=over 4

=item *

$configTxt: 

The configuration in plain text (either the file name or a reference to a string (in-memory file)

=back

I<Returns:>

Nothing

I<Exceptions:>

=over 4

Failed to open the configuration file

Failed to close the configuration file

Invalid syntax

=back

=cut

sub readTxt
{
    my ($self, $configTxt, $check) = @_;
    open my $handle, "<", $configTxt or errorExit("ERROR: Failed to open the configuration file $configTxt: $!");
    while(my $line = <$handle>)
    {
        chomp $line;
        $line =~ s/^\s+//; #remove leading spaces
        $line =~ s/\s+$//; #remove trailing spaces
        next unless $line;
        next if $line =~ /^#/;
        my @elementRefList = ($self);
        my $assignment = $line;
        if ($line =~ /^(\d+):\s*(\S.*)$/)
        {
            # set the variable for the given lanes
            $assignment = $2;
            @elementRefList = map($self->_getLane($_), split(//, $1));
        }
        # Check the the assignment is in the form "<NAME> <VALUE>"
        errorExit("ERROR: Invalid syntax: $assignment") unless $assignment =~ /^[A-Za-z][A-Za-z0-9_]*\s+\S+/;
        # Assign the variable
        foreach my $elementRef (@elementRefList)
        {
            $self->_assign($elementRef, $assignment, $check);
        }
    }
    close $handle or errorExit("ERROR: Failed to close the configuration file $configTxt: $!");
    # Provide a default EXPT_DIR if needed
    my $exptDir = $self->getVariable(name => 'EXPT_DIR');
    $self->setDefaultVariable('EXPT_DIR', '.') unless defined $exptDir and 0 < length($exptDir);
    # Make the EXPT_DIR absolute
    $self->setDefaultVariable('EXPT_DIR', File::Spec->rel2abs($self->getVariable(name => 'EXPT_DIR')));
}

sub readXml
{
    my ($self, $handler) = @_;
    my $xml = XMLin($handler, KeyAttr => [], ForceArray => ['Project', 'Reference', 'Sample', 'Barcode', 'Lane'], ForceContent => 1);
    foreach my $key (keys %$self)
    {
        $self->{$key} = $xml->{$key} if exists $xml->{$key};
    }
}

=head2 writeXml

Write the configuration as an XML file.

I<Parameters:>

=over 4

=item *

$configHandler

File handler to the config.xml

=back

I<Returns:>

Nothing.

I<Exceptions:>

=over 4

Failed to write the XML configuration

=back

=cut

sub writeXml
{
    my ($self, $file) = @_;
    open my $handler, '>', $file or errorExit("ERROR: Failed to open output file $file: $!");
    print $handler '<?xml version="1.0"?>'
    , "\n"
    , XMLout($self, RootName => 'RunParameters', KeyAttr => []) or errorExit("ERROR: Failed to write the XML configuration: $!");
    close $handler or errorExit("ERROR: Failed to close output file $file: $!");
}

=head2 toXml

Convert the configuration to the equivalent XML representation.

I<Parameters:>

None

I<Returns:>

The string containing the XML

I<Exceptions:>

=over 4

See 'writeXml'.

=back

=cut

sub toXml
{
    my ($self) = @_;
    my $xml;
    $self->writeXml(\$xml);
    return $xml;
}

=head2 getRunParameters

Load the run parameters from the EXPT_DIR DemultiplexedBustardConfig.xml.

I<Parameters:>

None

I<Returns:>

Nothing

I<Exceptions:>

=over 4

None.

=back

=cut

sub getRunParameters
{
    my $self = shift;
    errorExit("ERROR: The 'EXPT_DIR' must be specified") unless exists $self->{Defaults}->{EXPT_DIR}->{content};
    my $exptDir = $self->{Defaults}->{EXPT_DIR}->{content};
    errorExit("ERROR: The directory specified in 'EXPT_DIR' does not exist: $exptDir") unless -d $exptDir;
    my $configFile = File::Spec->catfile($exptDir, 'DemultiplexedBustardConfig.xml');
    my $xml = XMLin($configFile, 
                    KeyAttr => [], 
                    ForceArray => ['Lane', 'Reads', 'Tile', 'TileRange'],
                    ForceContent => 0) or errorExit("ERROR: Failed to load the XML config file $configFile");
    errorExit("ERROR: No 'Run' found in $configFile") unless exists $xml->{Run};
    my $run = $xml->{Run};
    errorExit("ERROR: No 'RunParameters' found in $configFile") unless exists $run->{RunParameters};
    my $runParameters = $run->{RunParameters};
    errorExit("ERROR: No 'Reads' element found in $configFile") unless exists $runParameters->{Reads};
    my $reads = $runParameters->{Reads};
    my @originalReads;
    foreach my $read (@$reads)
    {
        errorExit("ERROR: Missing 'Index' element for 'Reads' in $configFile") unless exists $read->{Index};
        errorExit("ERROR: Missing 'FirstCycle' element for 'Reads' in $configFile") unless exists $read->{FirstCycle};
        errorExit("ERROR: Missing 'LastCycle' element for 'Reads' in $configFile") unless exists $read->{LastCycle};
        my $index = $read->{Index};
        push @originalReads, $index;
        my $key = "ORIGINAL_READ_LENGTH$index";
        my $length = $read->{LastCycle} - $read->{FirstCycle} + 1;
        $self->{Defaults}->{$key}->{content} = $length unless exists $self->{Defaults}->{$key}->{content};
    }
    $self->{Defaults}->{ORIGINAL_READS}->{content} = join(' ', sort @originalReads) unless exists $self->{Defaults}->{ORIGINAL_READS}->{content};
    $self->{Defaults}->{FLOWCELL}->{content} = getFlowCellId($configFile) unless $self->{Defaults}->{FLOWCELL}->{content};

    my $runTileSelection = $run->{TileSelection};
    errorExit("ERROR: No 'TileSelection/Lane' element found in $configFile") unless exists $runTileSelection->{Lane};
    my $lanes = $runTileSelection->{Lane};
    foreach my $lane (@$lanes)
    {
        my @tiles = ();
        my $laneNumber = $lane->{Index};
        if (defined $lane->{Tile})
        {
            push @tiles, @{$lane->{Tile}};
        }
        if (defined $lane->{TileRange})
        {
            foreach my $tileRange (@{$lane->{TileRange}})
            {
                errorExit("ERROR: missing 'Min' in tile range for lane $laneNumber") unless exists $tileRange->{Min};
                errorExit("ERROR: missing 'Max' in tile range for lane $laneNumber") unless exists $tileRange->{Max};
                my $min = $tileRange->{Min};
                my $max = $tileRange->{Max};
                push @tiles, ($min..$max);
            }
        }
        if (@tiles)
        {
            my @sortedTiles = sort { $a <=> $b } @tiles;
            my $configLane = $self->_getLane($lane->{Index});
            $self->_setVariable('Tiles', join(' ', @sortedTiles), $configLane->{Defaults});
        }
    }
}

=head2 loadDefaults

Load the dictionary of default variables.

I<Parameters:>

None

I<Returns:>

Nothing

I<Exceptions:>

=over 4

None.

=back

=cut

sub loadDefaults
{
    my ($self, $handle) = @_;
    my $xml = XMLin($handle, SuppressEmpty => '');
    my $parameters = $xml->{DefaultRunParameters};
    my $targetElement = $self->{Defaults};
    foreach my $name (sort keys %$parameters)
    {
        my $value = $parameters->{$name};
        $self->_setVariable($name, $value, $targetElement);
    }
}


=head2 getVariable

Get the value of the given variables.

I<Parameters:>

(variable => <name>, lane => <index>, section => <Project|Reference|Sample|Barcode>, which => <section>)

I<Returns:>

The value of the variable for the required lane and section if any.

I<Exceptions:>

=over 4

None.

=back

=cut

sub getVariable
{
    my $self = shift;
    my %options = @_;
    errorExit("Missing variable name") unless exists $options{name};
    errorExit("Empty variable name") if $options{name} eq '';
    my $root = $self;
    $root = $self->_getLane($options{lane}) if exists $options{lane};
    my $selectedElement = $root->{Defaults};
    if (exists $options{section})
    {
        my $section = $options{section};
        errorExit("ERROR: Invalid section: $section") unless grep /^$section$/, ('Project', 'Reference', 'Sample', 'Barcode', 'Tile');
        $selectedElement = $self->_getElement($options{which}, $root, $section, 'name');
    }
    return undef unless exists $selectedElement->{$options{name}};
    my $variable = $selectedElement->{$options{name}};
    # empty variables don't have any content
    return '' unless %$variable;
    errorExit("Missing 'content' element") unless exists $variable->{content};
    return $variable->{content};
}

=head2 selectValue

Get the value of the given variables for the given dataset defined by (project, sample, lane, barcode).

The search is first at the 'Lanes' level then at the 'FlowCell' level. In both cases, the search order is
'Barcode', 'Sample', 'Project' and finally 'Default'.

I<Parameters:>

=over 4

=item *

$variable: The name of the variable

$project: the project for the dataset

$sample: the sample for the dataset

$lane: the lane for the dataset

$barcode: the barcode for the dataset

=back

I<Returns:>

The value of the variable for the specified dataset.

I<Exceptions:>

=over 4

None.

=back

=cut

sub selectValue
{
    my $self = shift;
    my ($variable, $project, $sample, $lane, $barcode, $reference) = @_;
    errorExit("ERROR: Missing variable name") unless defined $variable and $variable ne '';
    errorExit("ERROR: Missing project name") unless defined $project and $project ne '';
    errorExit("ERROR: Missing sample name") unless defined $sample and $sample ne '';
    errorExit("ERROR: Missing lane name") unless defined $lane and $lane ne '';
    errorExit("ERROR: Missing barcode name") unless defined $barcode and $barcode ne '';
    my %selection = (Barcode => $barcode, Sample => $sample, Project => $project, Reference => $reference);
    my @sections = qw(Barcode Sample Project);
    push @sections, 'Reference' if $reference;
    foreach my $root ($self->_getLane($lane), $self)
    {
        foreach my $section (@sections)
        {
            my $element = $self->_getElement($selection{$section}, $root, $section, 'name');
            next unless $element and exists $element->{$variable};
            return $element->{$variable}->{content} if exists $element->{$variable}->{content};
            # empty elements have no content
            return '';
        }
        next unless exists $root->{Defaults}->{$variable};
        return $root->{Defaults}->{$variable}->{content} if exists $root->{Defaults}->{$variable}->{content};
        # empty elements have no content
        return '';
    }
    errorExit("ERROR: Unknown configuration variable: $variable");
}

=head2 getDefaultVariables

Get the dictionary of default variables.

I<Parameters:>

None

I<Returns:>

The hash ( <VARIABLE> => <VALUE> ) with all the default variables

I<Exceptions:>

=over 4

None.

=back

=cut

sub getDefaultVariables
{
    my ($self) = @_;
    my %result;
    foreach my $name (sort keys %{$self->{Defaults}})
    {
        my $content = $self->{Defaults}->{$name};
        $result{$name} = undef;
        next unless $content and %$content;
        errorExit("ERROR: Missing 'content' element for default variable $name") unless $content and exists $content->{content};
        $result{$name} = $content->{content};
    }
    return %result;
}

sub setDefaultVariable
{
    my ($self, $name, $value, $check) = @_;
    $self->_setVariable($name, $value, $self->{Defaults}, $check);
}

=head2 validate

Validate the configuration

I<Parameters:>

None

I<Returns:>

Nothing

I<Exceptions:>

=over 4

Invalid configuration

=back

=cut

sub validate
{
    my ($self) = @_;
    foreach my $lane (@{$self->{Lanes}->{Lane}})
    {
        $self->_validate($lane);
    }
    $self->_validate($self);
}

=head2 validateDataset

Validate the configuration for a given dataset

I<Parameters:>

I<Parameters:>

=over 4

=item *

$project:

$sample:

$lane:

$barcode:

$reference:

=back

I<Returns:>

Nothing

I<Exceptions:>

=over 4

Invalid configuration

=back

=cut

sub validateDataset($$$$$)
{
    my $self = shift;
    my ($project, $sample, $lane, $barcode, $reference) = @_;
    my $analysis = $self->_validateExistence('ANALYSIS', $project, $sample, $lane, $barcode, $reference);
    my %validationMap = (none => \&_validateAnalysisNone,
                         eland_extended => \&_validateAnalysisElandExtended,
                         eland_pair => \&_validateAnalysisElandPair,
                         eland_rna => \&_validateAnalysisElandRna,
                         sequence => \&_validateAnalysisSequence);
    my @supportedAnalyses = keys %validationMap;
    errorExit("ERROR: unsupported analysis for $project-$sample-$lane-$barcode-$reference: $analysis: supported analyses are " .
              join(', ', @supportedAnalyses)) unless grep(/^$analysis$/, @supportedAnalyses);
    $validationMap{$analysis}->($self, $project, $sample, $lane, $barcode, $reference);
}

sub _validateExistence
{
    my $self = shift;
    my ($variable, $project, $sample, $lane, $barcode, $reference) = @_;
    errorExit("ERROR: no variable specified for $project-$sample-$lane-$barcode-$reference") unless $variable;
    my $value = $self->selectValue($variable, $project, $sample, $lane, $barcode, $reference);
    errorExit("ERROR: $variable not specified for $project-$sample-$lane-$barcode-$reference") unless $value;
    return $value;
}

sub _validatePath
{
    my $self = shift;
    my ($variable, $project, $sample, $lane, $barcode, $reference) = @_;
    my $value = $self->_validateExistence($variable, $project, $sample, $lane, $barcode, $reference);
    errorExit("ERROR: $variable for $project-$sample-$lane-$barcode-$reference: $value: not found") unless -e $value;
    return $value;
}

sub _validateDirectory
{
    my $self = shift;
    my ($variable, $project, $sample, $lane, $barcode, $reference) = @_;
    my $value = $self->_validatePath($variable, $project, $sample, $lane, $barcode, $reference);
    errorExit("ERROR: $variable for $project-$sample-$lane-$barcode-$reference: $value: not a directory") unless -d $value;
    return $value;
}

sub _validateAnalysisNone
{
    my $self = shift;
    my ($project, $sample, $lane, $barcode, $reference) = @_;
}
=com
sub _validateEland
{
    my $self = shift;
    my ($project, $sample, $lane, $barcode, $reference) = @_;
    for my $variable qw(ELAND_FASTQ_FILES_PER_PROCESS)
    {
        $self->_validateExistence($variable, $project, $sample, $lane, $barcode, $reference);
    }
    for my $directory qw(ELAND_GENOME)
    {
        $self->_validateDirectory($directory, $project, $sample, $lane, $barcode, $reference);
    }
}
=cut
sub _isBooleanTrue($)
{
    my $value = shift;
    return $value =~ /^y$|^yes$|^on$|^true$|^1$|^ok$/i;
}

sub _validateAnalysisElandExtended
{
    my $self = shift;
    my ($project, $sample, $lane, $barcode, $reference) = @_;
    $self->_validateEland($project, $sample, $lane, $barcode, $reference);

    my $multiRadsOverride = $self->selectValue('ELAND_EXTENDED_MULTI_READS', $project, $sample, $lane, $barcode, $reference);
    my $reads = $self->selectValue('READS', $project, $sample, $lane, $barcode, $reference);
    my @reads = split(/\s+/, $reads);
    errorExit("ERROR: eland_extended allows for one unmasked read only. Please adjust your USE_BASES for $project-$sample-$lane-$barcode-$reference") unless (1 == @reads or _isBooleanTrue($multiRadsOverride));
}


sub _validateAnalysisElandPair
{
    my $self = shift;
    my ($project, $sample, $lane, $barcode, $reference) = @_;
    $self->_validateEland($project, $sample, $lane, $barcode, $reference);

    my $reads = $self->selectValue('READS', $project, $sample, $lane, $barcode, $reference);
    my @reads = split(/\s+/, $reads);
    errorExit("ERROR: eland_pair requires two unmasked reads. Please adjust your USE_BASES for $project-$sample-$lane-$barcode-$reference") unless (2 == @reads);
}

sub _validateAnalysisElandRna
{
    my $self = shift;
    my ($project, $sample, $lane, $barcode, $reference) = @_;
    $self->_validateEland($project, $sample, $lane, $barcode, $reference);

    my $multiRadsOverride = $self->selectValue('ELAND_RNA_MULTI_READS', $project, $sample, $lane, $barcode, $reference);
    my $reads = $self->selectValue('READS', $project, $sample, $lane, $barcode, $reference);
    my @reads = split(/\s+/, $reads);
    errorExit("ERROR: eland_rna allows for one unmasked read only. Please adjust your USE_BASES for $project-$sample-$lane-$barcode-$reference") unless (1 == @reads or _isBooleanTrue($multiRadsOverride));

    $self->_validateExistence('ELAND_RNA_GENOME_ANNOTATION', $project, $sample, $lane, $barcode, $reference);
    my $value = $self->selectValue('ELAND_RNA_GENOME_ANNOTATION', $project, $sample, $lane, $barcode, $reference);

    my $geneMdGzSuffix = $self->selectValue('GENE_MD_GZ_SUFFIX', $project, $sample, $lane, $barcode, $reference);
    errorExit("ERROR: GENE_MD_GZ_SUFFIX must not be empty for eland_rna ANALYSIS") unless $geneMdGzSuffix;

    my $flatTxtGzSuffix = $self->selectValue('FLAT_TXT_GZ_SUFFIX', $project, $sample, $lane, $barcode, $reference);
    errorExit("ERROR: FLAT_TXT_GZ_SUFFIX must not be empty for eland_rna ANALYSIS") unless $flatTxtGzSuffix;

    errorExit("ERROR: ELAND_RNA_GENOME_ANNOTATION must end with $flatTxtGzSuffix or $geneMdGzSuffix") 
        unless $value =~ /$geneMdGzSuffix$|$flatTxtGzSuffix$/;
    $self->_validatePath('ELAND_RNA_GENOME_ANNOTATION', $project, $sample, $lane, $barcode, $reference);

    errorExit("ERROR: ELAND_RNA_GENE_MD_GROUP_LABEL must end be set") 
        unless $value !~ /$geneMdGzSuffix$/ || $self->selectValue('ELAND_RNA_GENE_MD_GROUP_LABEL', $project, $sample, $lane, $barcode, $reference);
}

sub _validateAnalysisSequence
{
    my $self = shift;
    my ($project, $sample, $lane, $barcode, $reference) = @_;
}

sub _validate
{
    my ($self, $root) = @_;
    foreach my $section ('Project', 'Reference', 'Sample', 'Barcode')
    {
        $self->_validateSection($root->{"${section}s"}->{$section}) if exists $root->{"${section}s"}->{$section};
    }
    $self->_validateSet($root->{Defaults});
}

sub _validateSection
{
    my ($self, $section) = @_;
    foreach my $set (@$section)
    {
        $self->_validateSet($set);
    }
}

sub _validateSet
{
    my ($self, $set) = @_;
    if (exists $set->{ANALYSIS})
    {
        my $content = $set->{ANALYSIS}->{content};
        errorExit("ERROR: Invalid ANALYSIS: $content") unless grep(/^$content$/, qw(none eland_extended eland_pair eland_rna sequence sequence_pair default));
    }
    if (exists $set->{EXPT_DIR})
    {
        my $content = $set->{EXPT_DIR}->{content};
        errorExit("ERROR: EXPT_DIR content does not exist") unless -e $content;
    }
    if (exists $set->{USE_BASES})
    {
        my @useBases = split(',', $set->{USE_BASES}->{content});
        my %originalReads = map {$_ => $self->{Defaults}->{"ORIGINAL_READ_LENGTH$_"}->{content}} split(/ /, $self->{Defaults}->{ORIGINAL_READS}->{content});
        push @useBases, (scalar(@useBases) ? $useBases[-1] : 'y*n') while scalar(@useBases) < scalar(keys %originalReads);
        my @reads;
        foreach my $index (sort keys %originalReads)
        {
            my $length = $originalReads{$index};
            my $mask = shift @useBases;
            my $expandedUseBases = expandUseBases($mask, $length);
            $set->{"USE_BASES$index"}->{content} = $expandedUseBases unless exists $set->{"USE_BASES$index"}->{content};
            my $maskedLength = $expandedUseBases =~ tr/(Y|y)/y/;
            $set->{"READ_LENGTH$index"}->{content} = $maskedLength unless exists $set->{"READ_LENGTH$index"}->{content};
            $set->{"ELAND_SEED_LENGTH$index"}->{content} = min($set->{"READ_LENGTH$index"}->{content}, 32)
                    unless (exists $set->{"ELAND_SEED_LENGTH$index"}->{content} && $set->{"ELAND_SEED_LENGTH$index"}->{content});
            push @reads, $index if $maskedLength;
        }
        $set->{READS}->{content} = join(' ', @reads) unless exists $set->{READS}->{content};
    }

    if (exists $set->{SAMTOOLS_GENOME} and $set->{SAMTOOLS_GENOME}->{content})
    {
        my $content = $set->{SAMTOOLS_GENOME}->{content};
        errorExit("ERROR: SAMTOOLS_GENOME ($content) file does not exist") unless -f $content;
        my ($volume, $dirs, $file) = File::Spec->splitpath($content);
        $set->{ELAND_GENOME}->{content} = File::Spec->canonpath(File::Spec->catpath($volume, $dirs));
        $set->{ELAND_GENOME_MASK}->{content} = $file;
        $set->{CHROM_NAME_SOURCE}->{content} = 'contigName';
    }
    elsif (exists $set->{ELAND_GENOME} and $set->{ELAND_GENOME}->{content})
    {
        $set->{ELAND_GENOME_MASK}->{content} = '*.fa' unless exists $set->{ELAND_GENOME_MASK}->{content};
        $set->{CHROM_NAME_SOURCE}->{content} = 'fileName';
        my $content = $set->{ELAND_GENOME}->{content};
        errorExit("ERROR: ELAND_GENOME ($content) directory does not exist") unless -d $content;
        errorExit("ERROR: ELAND_GENOME ($content) does not specify an absolute path") 
            unless File::Spec->file_name_is_absolute($content);
        my $elandGenomeGlobExpr = File::Spec->catfile($content, $set->{ELAND_GENOME_MASK}->{content});
        my @refFiles = glob($elandGenomeGlobExpr);
        errorExit("ERROR: ELAND_GENOME not a single file matches $elandGenomeGlobExpr") unless scalar(@refFiles);
    }
}

sub _assign
{
    my ($self, $elementRef, $assignment, $check) = @_;
    my $targetElement = $elementRef->{Defaults};
    foreach my $section (qw(Project Reference Sample Barcode))
    {
        my $ucSelection = uc($section);
        next unless $assignment =~ /^$ucSelection\s+(\S+)\s+(.*)$/;
        $targetElement = $self->_getElement($1, $elementRef, $section, 'name');
        $assignment = $2;
        last;
    }
    $assignment =~ /^(\S+)(\s+(\S.*)?)?$/;
    $self->_setVariable($1, $3, $targetElement, $check);
}

sub _setVariable
{
    my ($self, $name, $value, $root, $check) = @_;
    errorExit("ERROR: Undefined variable name") unless defined $name;
    if ($check)
    {
        errorExit("ERROR: Unknown variable: $name") unless exists $self->{Defaults}->{$name};
        # Some variables can be specified only for the whole flow cell
        my @flowCellVariables = qw(EXPT_DIR OUT_DIR);
        my $isFlowCellVariable = grep /^$name$/, @flowCellVariables;
        if ($root != $self->{Defaults} and $isFlowCellVariable)
        {
            errorExit("ERROR: Variable $name can be set only for the whole flow-cell");
        }
    }
    $root->{$name}->{content} = $value;
}

sub _getLane($$)
{
    my ($self, $index) = @_;
    my $laneRef = $self->_getElement($index, $self, 'Lane', 'index');
    $laneRef->{Defaults} = {} unless $laneRef->{Defaults};
    return $laneRef;
}

sub _getElement
{
    my ($self, $key, $root, $collectionName, $keyAttribute) = @_;
    errorExit("ERROR: Undefined collection name") unless defined $collectionName;
    errorExit("ERROR: Unknown collection name $collectionName") unless grep /^${collectionName}s$/, keys %$self;
    errorExit("ERROR: Undefined key attribute for $collectionName") unless defined $keyAttribute;
    $root->{"${collectionName}s"} = {$collectionName => []} unless $root->{"${collectionName}s"};
    my $collection = $root->{"${collectionName}s"}->{$collectionName};
    for my $elementRef (@$collection)
    {
        errorExit("ERROR: Missing key attribute $keyAttribute") unless defined $elementRef->{$keyAttribute};
        return $elementRef if $elementRef->{$keyAttribute} eq $key;
    }
    # Element not found. Create it
    my $elementRef = {$keyAttribute => $key};
    push @$collection, $elementRef;
    return $elementRef;
}

__END__

=pod

=head1 CONFIGURATION AND ENVIRONMENT

None

=head1 DEPENDENCIES

=over 4

=item Standard perl modules

strict, warnings, Exporter

=item External perl modules

XML::Simple

=item Casava perl modules

Casava::Common::Log

None

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Come Raczy

=cut
