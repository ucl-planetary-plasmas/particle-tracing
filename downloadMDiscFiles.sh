#!/bin/sh

#
# $Id: downloadMDiscFiles.sh,v 1.1 2017/11/13 11:02:24 patrick Exp $
#
# Copyright (c) 2017 Patrick Guio <patrick.guio@gmail.com>
# All Rights Reserved.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


# UCL model service site
UCL="http://astroweb.projects.phys.ucl.ac.uk/models/"

# -N Turn on time-stamping for wget
# -# Progress meter for curl
CMD="curl -#"

files="jup_mdisc_kh3e7_rmp90.mat sat_mdisc_kh2e6_rmp25.mat"
for f in $files; do
  if [ ! -e $f ] ; then
    echo "Downloading UCL Magnetodisc file $f"
		$CMD $UCL/mdisc/output/$f -o $f
	fi
done
