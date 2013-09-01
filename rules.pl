submit_rule(submit(CR, V)) :-
    base(CR, V),
    not(gerrit:commit_message_matches('\\\[WIP\\\]')),
    !.

submit_rule(submit(CR, V)) :-
    base(CR, V),
    not(gerrit:commit_message_matches('\\\[RFC\\\]')),
    !.

submit_rule(submit(CR, V)) :-
    base(CR, V).

base(CR, V) :-
    gerrit:max_with_block(-2, 2, 'Code-Review', CR),
    gerrit:max_with_block(-2, 2, 'Verified', V),
