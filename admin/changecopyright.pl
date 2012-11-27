#!/usr/bin/perl -w
#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2012, by the GROMACS development team, led by
# David van der Spoel, Berk Hess, Erik Lindahl, and including many
# others, as listed in the AUTHORS file in the top-level source
# directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org

# Use with
#
#  changecopyright.pl your_file_name_here
#
# Files with C-style comments and CMake files are more-or-less supported

use strict;
use Regexp::Common qw /comment RE_comment_C/;
use File::Slurp;

my $file_name = shift @ARGV;
my $gromacs_copyright_notice = <<END;
/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
END
chop $gromacs_copyright_notice;

my %filetypes = (
    [ 'extension'   => qr/.*\.(?:[chyl]|[ch]pp)|cuh?$|/,
      'commentchar' => ''
    ]
    ,
    [ 'extension'   => qr/.*\.f$/,
      'commentchar' => 'C'
    ]
    ,
    [ 'extension'   => qr/.*\.(?:py|pl|sh|)$/,
      'commentchar' => '#'
    ]
    ,
    [ 'extension'   => qr/CMakeLists.txt/,
      'commentchar' => '#'
    ]
    ,
    [ 'extension'   => qr/.pre/,
      'commentchar' => ''
    ]
    ,
    [ 'extension'   => qr//,
      'commentchar' => ''
    ]
    );

# since the new copyright notice is not included here, we magically
# can't duplicate the 4.6 copyright notice by re-applying the script
my $multi_line_authors = qr/ \* BIOSON Research Institute, Dept. of Biophysical Chemistry\n \* University of Groningen, The Netherlands\n| \* check out http:\/\/www.gromacs.org for more information.\n| \* Erik Lindahl, David van der Spoel, University of Groningen.\n| \* David van der Spoel, Erik Lindahl, University of Groningen.\n|\* David van der Spoel, Erik Lindahl, University of Groningen.\n| \* David van der Spoel, Erik Lindahl, Berk Hess, University of Groningen.\n/;

my $year = qr/[0-9]{4}/;
my $year_range = qr/$year(?:-$year)?/;
my $noncopyright_line = qr/ \* (?![Cc]opyright)[^\n]*\n/;
my $empty_line = qr/ ?\*\n/;
my $text = qr/[^\n]/;
my $comma_and_more_lines = qr/,\n${multi_line_authors}/;
my $nothing_and_possible_multi_line_authors = qr/ *\n${multi_line_authors}?/;
my $period_and_end_of_line = qr/\.\s*?\n/;
my $copyright_notice = qr/(?:This file is part of Gromacs        |Gromacs 4.0                         )?[Cc]opyright.*\s(${year_range}(?:,${year_range})*)(?:${text}+${period_and_end_of_line}|${text}+${comma_and_more_lines}|${nothing_and_possible_multi_line_authors}|${text}+\n(?=${empty_line}))/;

my $slurped_text = read_file($file_name);
my $result = "";
my $remaining = $slurped_text;
my $found_copyright_notice = 0;

while ($remaining =~ /$RE{comment}{C}/) {
    # These variables are ugly, but so is the purpose of this script
    $result .= ${^PREMATCH};
    my $comment = ${^MATCH};
    $remaining = ${^POSTMATCH};

    # Accept this comment unmodified, and later prepend the GROMACS
    # one because found_copyright_notice is still 0
    if ($comment =~ /This source code file is part of thread_mpi|Aladdin Enterprises/) {
        $result .= $comment;
        next;
    }
    if ($comment =~ $copyright_notice && 0 == $found_copyright_notice) {
        my $this_copyright_notice = $gromacs_copyright_notice;
        my @old_copyrights = ();
        my $remaining_comment;
        do {
            push @old_copyrights, ${^MATCH};
            $remaining_comment = ${^POSTMATCH};
        } while ($remaining_comment =~ /${copyright_notice}/m);

        # Prepend old copyrights to the new GROMACS copyright notice
        map {
            (my $old_copyright = $_) =~ s/,$/./;
            # Remove temporary open-ended copyright assertions
            unless ($old_copyright =~ /Copyright \(c\) 2012- */) {
                $this_copyright_notice =~ s/^( \* Copyright)/ * ${old_copyright}$1/m;
            }
        } reverse(@old_copyrights);

        $result .= $this_copyright_notice;
    } else {
        $result .= $comment;
    }
    if ($comment =~ /[Cc]opyright/) {
        # Make sure we don't duplicate copyright statements when the first
        # one was made in 2012, but can still insert that statement when
        # there has never been one before.
        $found_copyright_notice = 1;
    }
}
$result .= $remaining;

# Deal with CMake files first
if ($file_name =~ /CMakeLists.txt$|.cmake$/) {
    if ( $result !~ /[Cc]opyright.*by the GROMACS development team/) {
        (my $this_copyright_notice = $gromacs_copyright_notice)  =~ s/^[\/ ]\*\/?/#/gm;
        $result = "${this_copyright_notice}\n$result";
    }
} elsif (0 == $found_copyright_notice) {
    # Prepare output for files with C-style comments
    $result = "${gromacs_copyright_notice}\n$result";
}

# Write output
open FILE, ">$file_name" or die "Couldn't open file $file_name for writing";
print FILE $result;
close FILE;

# Review the changes before they get added! Some files need manual tweaking
# because they have licensing external to GROMACS.
system("git add --interactive --patch $file_name");
