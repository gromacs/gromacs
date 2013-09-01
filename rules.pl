submit_rule(submit(CR, V, RFC)) :-
    base(CR, V),
    gerrit:commit_message_matches('\\\[WIP\\\]|\\\[RFC\\\]'),
    !,
    RFC = label('RFC-or-WIP', ok(_)).

submit_rule(submit(CR, V, RFC)) :-
    base(CR, V),
    RFC = label('RFC-or-WIP', need(_)).

base(CR, V) :-
    gerrit:max_with_block(-2, 2, 'Code-Review', CR),
    gerrit:max_with_block(-2, 2, 'Verified', V).
