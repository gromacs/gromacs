Logging {#page_logging}
=======

Currently, mdrun is using a combination of direct C-style I/O into `fplog` and
`stderr`, and the facilities described here.  However, more and more should get
moved to this interface in the future.

The parts that make up the logging system are shown below.

\dot
    digraph logging_overview {
        builder [label="LoggerBuilder", URL="\ref gmx::LoggerBuilder"]
        owner [label="LoggerOwner", URL="\ref gmx::LoggerOwner"]
        logger [label="MDLogger", URL="\ref gmx::MDLogger"]
        target [label="ILogTarget", URL="\ref gmx::ILogTarget"]
        user [label="using code"]

        builder -> owner [label="builds"]
        owner -> logger
        owner -> target [label="owns"]
        logger -> target [label="references"]
        user -> builder [label="set logging targets"]
        user -> logger [label="write with\nGMX_LOG()"]
    }
\enddot

To initialize the logging system, the using code creates an instance of
gmx::LoggerBuilder, and sets the desired logging targets with provided methods.
Once all targets have been initialized, the code calls
gmx::LoggerBuilder::build() and gets a gmx::LoggerOwner, which is responsible
of managing the memory allocated for the logger.

To log information, the using code uses an gmx::MDLogger returned by
gmx::LoggerOwner::logger() with the ::GMX_LOG macro.  Code that writes to the
log only needs to know of this class (and helper classes used to implement the
macro), which is a relatively simple container for references to the logging
targets.  If there is no log target that would consume the information written
with ::GMX_LOG, the whole statement evaluates to a conditional that reads the
log target from a member variable and compares it against `nullptr`.  All the
code that formats the output is skipped in this case.

Currently the implementation is geared to making ::GMX_LOG behavior stable, and
to be relatively extensible.  However, using any other approach than ::GMX_LOG
for writing to the log should first think about how the API could be best
organized for that.

All information written to the log is composed of _log entries_.  Each
::GMX_LOG statement writes a single log entry, meaning that newlines are
automatically added.

The logging methods are not thread-safe, so it is the responsibility of the
calling code to only use them from a single thread or otherwise synchronize
access.
