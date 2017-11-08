#!/usr/bin/perl
#Script to test authenticated resources via REST interface.
#Written by Keith Jolley
#Modified by Adam Koziol
#Original Copyright (c) 2015, University of Oxford
#E-mail: keith.jolley@zoo.ox.ac.uk
#
#This file is part of Bacterial Isolate Genome Sequence Database (BIGSdb).
#
#BIGSdb is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#BIGSdb is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#Usage: rest_auth.pl <route>
use strict;
use warnings;
use 5.010;
use Net::OAuth 0.20;
$Net::OAuth::PROTOCOL_VERSION = Net::OAuth::PROTOCOL_VERSION_1_0A;
use HTTP::Request::Common;
use LWP::UserAgent;
use JSON qw(encode_json decode_json);
use Data::Random qw(rand_chars);
use Data::Dumper;
use Config::Tiny;
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;
use POSIX;
use MIME::Base64;
use Time::Piece;

# This start time will be used in calculating the total time of the run
my $start_time = time;

#Modify the following values to your local system.
use constant TEST_REST_URL   => 'http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef';
use constant TEST_WEB_URL    => 'http://pubmlst.org/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_rmlst_seqdef';
###
my %opts;
GetOptions(
	'a|arguments=s'     => \$opts{'a'},
	'f|file=s'          => \$opts{'f'},
	'i|isolates_file=s' => \$opts{'i'},
	'm|method=s'        => \$opts{'m'},
	'p|profiles_file=s' => \$opts{'p'},
	'r|route=s'         => \$opts{'r'},
	's|sequence_file=s' => \$opts{'s'},
	'h|help'            => \$opts{'h'},
	'prompt'            => \$opts{'prompt'}
) or die("Error in command line arguments\n");
if ( $opts{'h'} ) {
	show_help();
	exit;
}

### Edits - open the secrets file containing the consumer key and consumer secret
my $secretfile = $opts{'a'};

# Open the file
open my $fh, '<', $secretfile or die "Cannot find the secret.txt file required for authorization. Please
ensure that this file exists, and that the supplied consumer key is on the first line, and the consumer secret is on 
the second line. Contact keith.jolley\@zoo.ox.ac.uk for an account, and the necessary keys.\n";
# Read the file to an array
my @data = <$fh>;
# Close the file
close $fh;

chomp (@data);
my $consumerkey = $data[0];
my $consumersecret = $data[1];
###


my $ua      = LWP::UserAgent->new;
my @methods = qw(GET POST PUT DELETE);
my %methods = map { $_ => 1 } @methods;
$opts{'m'} //= 'GET';
local $" = q(, );
die "Invalid method - allowed values are: @methods.\n" if !$methods{ $opts{'m'} };
main();

# Return the run time
my $end_time = time;
my $total_time = $end_time - $start_time;
my $legible_time = sprintf("%.1d", $total_time);

print "-------------------------------------------\n";
print "The total run time was $legible_time seconds.\n\n";

#Access and session tokens are stored within current directory.
#If a session token is missing or expired, a new one will be requested using the access token.
#If an access token is missing or expired, a new one will be requested.
sub main {
	my $route = $opts{'r'} // q();
	my ( $session_token, $session_secret ) = _retrieve_token('session_token');
	if ( !$session_token || !$session_secret ) {
		my $session_response = _get_session_token();
		( $session_token, $session_secret ) = ( $session_response->token, $session_response->token_secret );
	}
	_get_route( $route, $session_token, $session_secret );
	return;
}

sub _get_request_token {
	my $request = Net::OAuth->request('request token')->new(
		consumer_key     => $consumerkey,
		consumer_secret  => $consumersecret,
		request_url      => TEST_REST_URL . '/oauth/get_request_token',
		request_method   => 'GET',
		signature_method => 'HMAC-SHA1',
		timestamp        => time,
		nonce            => join( '', rand_chars( size => 16, set => 'alphanumeric' ) ),
		callback         => 'oob'
	);
	$request->sign;

	#say $request->signature_base_string;
	die "COULDN'T VERIFY! Check OAuth parameters.\n" unless $request->verify;
	say 'Getting request token...';
	_prompt();
	my $res = $ua->request( GET $request->to_url, Content_Type => 'application/json', );
	my $decoded_json = decode_json( $res->content );
	my $request_response;
	if ( $res->is_success ) {
		say 'Success:';
		$request_response = Net::OAuth->response('request token')->from_hash($decoded_json);
		say 'Request Token:        ', $request_response->token;
		say 'Request Token Secret: ', $request_response->token_secret;
		_write_token( 'request_token', $request_response->token, $request_response->token_secret );
		return $request_response;
	} else {
		say 'Failed:';
		exit;
	}
}

sub _get_access_token {
	my ( $request_token, $request_secret ) = @_;
	unlink 'access_token';
	if ( !$request_token || $request_secret ) {
		( $request_token, $request_secret ) = _retrieve_token('request_token');
		if ( !$request_token || !$request_secret ) {
			my $session_response = _get_request_token();
			( $request_token, $request_secret ) = ( $session_response->token, $session_response->token_secret );
		}
	}
	say "\nNow log in at\n"
	  . TEST_WEB_URL
	  . "&page=authorizeClient&oauth_token=$request_token"
	  . "\nto obtain a verification code.";
	print "\nPlease enter verification code:  ";
	my $verifier = <>;
	chomp $verifier;
	my $request = Net::OAuth->request('access token')->new(
		consumer_key     => $consumerkey,
		consumer_secret  => $consumersecret,
		token            => $request_token,
		token_secret     => $request_secret,
		verifier         => $verifier,
		request_url      => TEST_REST_URL . '/oauth/get_access_token',
		request_method   => 'GET',
		signature_method => 'HMAC-SHA1',
		timestamp        => time,
		nonce            => join( '', rand_chars( size => 16, set => 'alphanumeric' ) ),
	);
	$request->sign;
	die "COULDN'T VERIFY! Check OAuth parameters.\n" unless $request->verify;
	say "\nGetting access token...";
	unlink 'request_token';    #Request tokens can only be redeemed once
	my $res = $ua->request( GET $request->to_url, Content_Type => 'application/json' );
	my $decoded_json = decode_json( $res->content );

	if ( $res->is_success ) {
		say 'Success:';
		my $access_response = Net::OAuth->response('access token')->from_hash($decoded_json);
		say 'Access Token:        ', $access_response->token;
		say 'Access Token Secret: ', $access_response->token_secret;
		say "\nThis access token will not expire but may be revoked";
		say 'by the user or the service provider. It may be used to';
		say 'obtain temporary session tokens.';
		_write_token( 'access_token', $access_response->token, $access_response->token_secret );
		return $access_response;
	} else {
		say 'Failed:';
		return;
	}
}

sub _get_session_token {
	my ( $access_token, $access_secret ) = @_;
	unlink 'session_token';
	if ( !$access_token || $access_secret ) {
		( $access_token, $access_secret ) = _retrieve_token('access_token');
		if ( !$access_token || !$access_secret ) {
			my $session_response = _get_access_token();
			( $access_token, $access_secret ) = ( $session_response->token, $session_response->token_secret );
		}
	}
	say "\nNow requesting session token using access token...";
	_prompt();
	my $request = Net::OAuth->request('protected resource')->new(
		consumer_key     => $consumerkey,
		consumer_secret  => $consumersecret,
		token            => $access_token,
		token_secret     => $access_secret,
		request_url      => TEST_REST_URL . '/oauth/get_session_token',
		request_method   => 'GET',
		signature_method => 'HMAC-SHA1',
		timestamp        => time,
		nonce            => join( '', rand_chars( size => 16, set => 'alphanumeric' ) ),
	);
	$request->sign;

	#say $request->signature_base_string;
	die "COULDN'T VERIFY! Check OAuth parameters.\n" unless $request->verify;
	say "\nGetting session token...";
	my $res = $ua->request( GET $request->to_url, Content_Type => 'application/json' );
	my $decoded_json = decode_json( $res->content );
	if ( $res->is_success ) {
		say 'Success:';
		my $session_response = Net::OAuth->response('access token')->from_hash($decoded_json);
		say 'Session Token:        ', $session_response->token;
		say 'Session Token Secret: ', $session_response->token_secret;
		say "\nThis session token will expire in 12 hours (default).";
		say 'It should be used with the secret to sign any requests';
		say 'to the API.';
		_write_token( 'session_token', $session_response->token, $session_response->token_secret );
		return $session_response;
	} else {
		say 'Failed:';
		if ( $res->{'_content'} =~ /401/ ) {
			say 'Invalid access token, requesting new one...';
			my $access_response = _get_access_token();
			if ($access_response) {
				( $access_token, $access_secret ) = ( $access_response->token, $access_response->token_secret );
			}
			return _get_session_token( $access_token, $access_secret );
		} else {
			return;
		}
	}
}

sub _get_route {
	my ( $route, $session_token, $session_secret ) = @_;
	$route //= q();
	if ($route) {
		$route = "/$route" if $route !~ /^\//x;
	}
	my $extra_params = {};
	if ( $opts{'s'} ) {
		die "Sequence file $opts{'s'} does not exist.\n" if !-e $opts{'s'};
		my $seqs = _slurp( $opts{'s'} );
		$extra_params->{'sequences'} = $$seqs;
	}
	if ( $opts{'p'} ) {
		die "Profiles file $opts{'p'} does not exist.\n" if !-e $opts{'p'};
		my $profiles = _slurp( $opts{'p'} );
		$extra_params->{'profiles'} = $$profiles;
	}
	if ( $opts{'i'} ) {
		die "Isolates file $opts{'i'} does not exist.\n" if !-e $opts{'i'};
		my $isolates = _slurp( $opts{'i'} );
		$extra_params->{'isolates'} = $$isolates;
	}
	if ( $opts{'f'} ) {
		die "File $opts{'f'} does not exist.\n" if !-e $opts{'f'};
		open( my $fh, '<', $opts{'f'} ) || die "Cannot open file $opts{'f'} for reading.\n";
		binmode $fh;
		my $contents = do { local $/ = undef; <$fh> };
		close $fh;
		$extra_params->{'upload'} = encode_base64($contents);
	}
	if ( $opts{'a'} ) {
		my @p = split /&/x, $opts{'a'};
		for my $pa (@p) {
			my @ps = split /=/x, $pa;
			$extra_params->{ $ps[0] } = $ps[1];
		}
	}
	my $url = TEST_REST_URL . "$route";
	say "\nAccessing authenticated resource ($url)...";
	my $request = Net::OAuth->request('protected resource')->new(
		consumer_key     => $consumerkey,
		consumer_secret  => $consumersecret,
		token            => $session_token,
		token_secret     => $session_secret,
		# Edit
		request_url      => $url . "/loci",
		request_method   => $opts{'m'},
		signature_method => 'HMAC-SHA1',
		timestamp        => time,
		nonce            => join( '', rand_chars( size => 16, set => 'alphanumeric' ) ),
		extra_params     => $extra_params
	);
	$request->sign;
	die "COULDN'T VERIFY! Check OAuth parameters.\n" unless $request->verify;
	my $request_params = $request->all_message_params;
	my %all_params     = %$extra_params;
	foreach my $param (@$request_params) {
		$all_params{"oauth_$param"} = $request->$param;
	}
	my $method = lc( $opts{'m'} );
	my $res =
	    $opts{'m'} eq 'POST'
	  ? $ua->post( $url, Content => encode_json \%all_params )
	  : $ua->$method( $request->to_url );
	my $decoded_json;
	eval { $decoded_json = decode_json( $res->content ) };
	if ($@) {
		say $res->content;
		return;
	}
	
	
	if ( ( $decoded_json->{'message'} // q() ) =~ /Client\ is\ unauthorized/x ) {
		say 'Access denied - client is unauthorized.';
		return;
	}
	if ( ( $decoded_json->{'status'} // q() ) eq '401' ) {
		say 'Invalid session token, requesting new one...';
		my $session_response = _get_session_token();
		if ($session_response) {
			( $session_token, $session_secret ) = ( $session_response->token, $session_response->token_secret );
		}
		_get_route( $route, $session_token, $session_secret );
	}
	
	# Edits	
	_get_route ( $route, $session_token, $session_secret ) if ref $decoded_json ne 'HASH';
	say Dumper ($decoded_json);
	# Get list of loci from the request
	my @loci;
	foreach my $post ($decoded_json) {
		foreach my $loci ($post->{'loci'}) {
			# Set @loci to the list of loci discovered
			@loci = @$loci;
		}
	}
	# Download profile
	my $path = getcwd;
#	mkdir ("$path/rMLST");
	downloadprofile($session_token, $session_secret, $opts{'m'}, $extra_params, $url, $route, $path);
	
	# Download the alleles
	use threads;
	my @threads;
	# Determine the number of threads present in the system
	my @cpus = `awk '/^processor/ { N++} END { print N }' /proc/cpuinfo`;
	chomp @cpus;
	# Multi-threaded approach to format the database
	while (scalar(@threads) < @cpus) {
		foreach my $locus (@loci) {
			my $t = threads->new(\&downloadallele, $locus, $session_token, $session_secret, $opts{'m'}, $extra_params, $url, $route, $path);
			# Output data to @threads. Right now, all this does is increase the size for the scalar(@threads) portion of the while loop
		push(@threads,$t);
		}
	}
	
	# This loop ensures that each thread is complete before terminating
	foreach (@threads) {
		# The join command is responsible for ensuring that all threads are complete
		my $num = $_->join;
	}
	
	return;
}


##########################

sub downloadprofile {
	# Download the profile
	my ($session_token, $session_secret, $opts, $extra_params, $url, $route, $path) = @_;
	# Define the file name for the allele file
	my $profile = "$path/profile.txt";
	# Find the size of the profile file
	my $profilesize = (stat $profile)[7];
	# Set the size variable to 0 if stat could not determine the size
	unless ($profilesize){
		$profilesize = 0
	}
	# Make the request to download the profile if the size of the profile is under 100 bytes
	if ($profilesize < 100) {
		say "Downloading profile";
		my $url = TEST_REST_URL . "$route";
		my $request = Net::OAuth->request('protected resource')->new(
			consumer_key     => $consumerkey,
			consumer_secret  => $consumersecret,
			token            => $session_token,
			token_secret     => $session_secret,
			request_url      => $url . "/schemes/1/profiles_csv",
			request_method   => $opts{'m'},
			signature_method => 'HMAC-SHA1',
			timestamp        => time,
			nonce            => join( '', rand_chars( size => 16, set => 'alphanumeric' ) ),
			extra_params     => $extra_params
		);
		$request->sign;
		#say $request->signature_base_string;
		die "COULDN'T VERIFY! Check OAuth parameters.\n" unless $request->verify;
		my $request_params = $request->all_message_params;
		my %all_params     = %$extra_params;
		foreach my $param (@$request_params) {
			$all_params{"oauth_$param"} = $request->$param;
		}
		my $method = lc( $opts{'m'} );
		my $res =
		    $opts{'m'} eq 'POST'
		  ? $ua->post( $url, Content => encode_json \%all_params )
		  : $ua->$method( $request->to_url );
		my $decoded_json;
		eval { $decoded_json = decode_json( $res->content ) };
		
		
#		return if ref $decoded_json ne 'HASH';
		if ( ( $decoded_json->{'message'} // q() ) =~ /Client\ is\ unauthorized/x ) {
			say 'Access denied - client is unauthorized.';
			return;
		}
		if ( ( $decoded_json->{'status'} // q() ) eq '401' ) {
			say 'Invalid session token, requesting new one...';
			my $session_response = _get_session_token();
			if ($session_response) {
				( $session_token, $session_secret ) = ( $session_response->token, $session_response->token_secret );
			}
			_get_route( $route, $session_token, $session_secret );
		}
		if ($@) {
			unless ($res->content =~ /read timeout/) {
				open(PROFILE, ">", $profile) or die $!;
				say PROFILE $res->content;
				close(PROFILE);
			}
			else{
				say "Profile did not download. Trying again.";
				downloadprofile($session_token, $session_secret, $opts, $extra_params, $url, $route, $path);
			}	
		
		} else {
			print "Profile did not download. Trying again.";
			downloadprofile($session_token, $session_secret, $opts, $extra_params, $url, $route, $path);
		}
	}
}


##############
	
sub downloadallele {
	my ($locus, $session_token, $session_secret, $opts, $extra_params, $url, $route, $path) = @_;
	# Split the name of the gene from the URL:
	my @split = split /\//, $locus;
	# Define the file name for the allele file
	my $fasta = "$path/$split[-1].tfa";
	# Use the stat method to determine the size of the fasta file	
	my $fastasize = (stat $fasta)[7];
	# Set $fastasize to 0 if there was an issue using stat
	unless ($fastasize){
		$fastasize = 0
	}
	# Only perform the request if the file size is less than 100 bytes
	if ($fastasize < 100) {
		my $request = Net::OAuth->request('protected resource')->new(
			consumer_key     => $consumerkey,
			consumer_secret  => $consumersecret,
			token            => $session_token,
			token_secret     => $session_secret,
			request_url      => $locus . "/alleles_fasta",
			request_method   => $opts{'m'},
			signature_method => 'HMAC-SHA1',
			timestamp        => time,
			nonce            => join( '', rand_chars( size => 16, set => 'alphanumeric' ) ),
			extra_params     => $extra_params
		);
		$request->sign;
		die "COULDN'T VERIFY! Check OAuth parameters.\n" unless $request->verify;
		my $request_params = $request->all_message_params;
		my %all_params     = %$extra_params;
		foreach my $param (@$request_params) {
			$all_params{"oauth_$param"} = $request->$param;
		}
		my $method = lc( $opts{'m'} );
		my $res =
		    $opts{'m'} eq 'POST'
		  ? $ua->post( $url, Content => encode_json \%all_params )
		  : $ua->$method( $request->to_url );
		my $decoded_json;
		eval { $decoded_json = decode_json( $res->content ) };
		# If data were received
		if ($@) {
			# If the request was successful
			unless ($res->content =~ /read timeout/) {
				say "Saving $split[-1] allele file";
				# Save the data to file
				open(ALLELEFILE, ">", $fasta) or die $!;
				say ALLELEFILE $res->content;
				close(ALLELEFILE);
			}
			# If there was a time out trying to retrieve the requested resource
			else{
				say "$fasta did not download. Trying again";
				downloadallele($locus, $session_token, $session_secret, $opts, $extra_params, $url, $route, $path);
			}	
		
		} else {
			say "$fasta did not download. Trying again";
			downloadallele($locus, $session_token, $session_secret, $opts, $extra_params, $url, $route, $path);
		}
	}
}


sub _retrieve_token {
	my ($token_name) = @_;
	return if !-e $token_name;
	my $config = Config::Tiny->new();
	$config = Config::Tiny->read($token_name);
	return ( $config->{_}->{'token'}, $config->{_}->{'secret'} );
}

sub _write_token {
	my ( $token_name, $token, $secret ) = @_;
	my $config = Config::Tiny->new();
	$config->{_}->{'token'}  = $token;
	$config->{_}->{'secret'} = $secret;
	$config->write($token_name);
	return;
}

sub _prompt {
	return if !$opts{'prompt'};
	say 'Press any key to continue...';
	<STDIN>;
	return;
}

sub _slurp {
	my ($file_path) = @_;
	open( my $fh, '<:encoding(utf8)', $file_path )
	  || throw BIGSdb::CannotOpenFileException("Can't open $file_path for reading");
	my $contents = do { local $/ = undef; <$fh> };
	return \$contents;
}

sub show_help {
	my $termios = POSIX::Termios->new;
	$termios->getattr;
	my $ospeed = $termios->getospeed;
	my $t = Tgetent Term::Cap { TERM => undef, OSPEED => $ospeed };
	my ( $norm, $bold, $under ) = map { $t->Tputs( $_, 1 ) } qw/me md us/;
	say << "HELP";
${bold}NAME$norm
    ${bold}rest_auth.pl$norm - Test authenticated resources via REST interface

${bold}SYNOPSIS$norm
    ${bold}rest_auth.pl --route$norm ${under}ROUTE$norm [${under}options$norm]

${bold}OPTIONS$norm

${bold}-a, --arguments$norm ${under}DATA$norm
    Data to put in during POST e.g. 
    'type=alleles&software=Enterobase'.

${bold}-f, --file$norm ${under}FILENAME$norm
    Name of file to upload.

${bold}-h, --help$norm
    This help page.
    
${bold}-i, --isolates_file$norm ${under}FILE$norm
    Relative path of tab-delimited file of isolate data to upload.
   
${bold}-m, --method$norm ${under}GET|PUT|POST|DELETE$norm
    Set HTTP method (default GET).
    
${bold}-p, --profiles_file$norm ${under}FILE$norm
    Relative path of tab-delimited file of allelic profiles to upload.
    
${bold}--prompt$norm
    Prompt before connection requests (used for demonstration purposes).

${bold}-r, --route$norm ${under}ROUTE$norm
    Relative path of route, e.g. 'submissions'.
    
${bold}-s, --sequence_file$norm ${under}FILE$norm
    Relative path of FASTA or single sequence file to upload.
HELP
	return;
}
