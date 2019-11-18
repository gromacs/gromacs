#!/usr/bin/perl -w

# in: input filename
$in = shift || die("Please specify input filename");
# If your exchange was every N ps and you saved every M ps you can make for
# the missing frames by setting extra to (N/M - 1). If N/M is not integer,
# you're out of luck and you will not be able to demux your trajectories at all.
$extra = shift || 0;
$ndx  = "replica_index.xvg";
$temp = "replica_temp.xvg";

@comm = ("-----------------------------------------------------------------",
	 "Going to read a file containing the exchange information from",
	 "your mdrun log file ($in).", 
	 "This will produce a file ($ndx) suitable for",
	 "demultiplexing your trajectories using trjcat,",
	 "as well as a replica temperature file ($temp).",
	 "Each entry in the log file will be copied $extra times.",
	 "-----------------------------------------------------------------");
for($c=0; ($c<=$#comm); $c++) {
    printf("$comm[$c]\n");
}

# Open input and output files
open (IN_FILE,"$in") || die ("Cannot open input file $in");
open (NDX,">$ndx") || die("Opening $ndx for writing");
open (TEMP,">$temp") || die("Opening $temp for writing");


sub pr_order {
    my $t     = shift;
    my $nrepl = shift;
    printf(NDX "%-20.2f",$t);
    for(my $k=0; ($k<$nrepl); $k++) {
	my $oo = shift;
	printf(NDX "  %3d",$oo);
    }
    printf(NDX "\n");
}

sub pr_revorder {
    my $t     = shift;
    my $nrepl = shift;
    printf(TEMP "%-20.2f",$t);
    for(my $k=0; ($k<$nrepl); $k++) {
	my $oo = shift;
	printf(TEMP "  %3d",$oo);
    }
    printf(TEMP "\n");
}

$nrepl = 0;
$init  = 0;
$tstep = 0;
$nline = 0;
$tinit = 0;
while ($line = <IN_FILE>) {
    chomp($line);
    
    if (index($line,"init_t") >= 0) {
	@log_line = split (' ',$line);
	$tinit = $log_line[2];
    }
    if (index($line,"Repl") == 0) {
	@log_line = split (' ',$line);
	if (index($line,"There") >= 0) {
	    $nrepl = $log_line[3];
	}
	elsif (index($line,"time") >= 0) {
	    $tstep = $log_line[6];
	}
	elsif ((index($line,"Repl ex") == 0) && ($nrepl == 0)) {
            # Determine number of replicas from the exchange information
	    printf("%s\n%s\n",
		   "WARNING: I did not find a statement about number of replicas",
		   "I will try to determine it from the exchange information.");
	    for($k=2; ($k<=$#log_line); $k++) {
		if ($log_line[$k] ne "x") {
		    $nrepl++;
		}
	    }
	}
	if (($init == 0) && ($nrepl > 0)) {
	    printf("There are $nrepl replicas.\n");

	    @order = ();
            @revorder = ();
	    for($k=0; ($k<$nrepl); $k++) {
		$order[$k] = $k;
                $revorder[$k] = $k;
	    }
	    for($ee=0; ($ee<=$extra); $ee++) {
		pr_order($tinit+$ee,$nrepl,@order);
		pr_revorder($tinit+$ee,$nrepl,@revorder);
		$nline++;
	    }
	    $init = 1;
	}

	if (index($line,"Repl ex") == 0) {
	    $k = 0;
	    for($m=3; ($m<$#log_line); $m++) {
		if ($log_line[$m] eq "x") {
		    $revorder[$order[$k]] = $k+1;
		    $revorder[$order[$k+1]] = $k;
		    $tmp = $order[$k];
		    $order[$k] = $order[$k+1];
		    $order[$k+1] = $tmp;
#	    printf ("Swapping %d and %d on line %d\n",$k,$k+1,$line_number); 
		}
		else {
		    $k++;
		}
	    }
	    for($ee=0; ($ee<=$extra); $ee++) {
		pr_order($tstep+$ee,$nrepl,@order);
		pr_revorder($tstep+$ee,$nrepl,@revorder);
		$nline++;
	    }
	}
    }
}
close IN_FILE;
close NDX;
close TEMP;

printf ("Finished writing $ndx and $temp with %d lines\n",$nline);
