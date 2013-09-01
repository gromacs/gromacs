submit_rule(submit(CR, V, RFC)) :-
    base(CR, V),
    gerrit:commit_message_matches('(\\\[WIP\\\]|\\\[RFC\\\])'),
    !,
    RFC = label('Not-RFC-or-WIP', need(_)).

submit_rule(submit(CR, V, RFC)) :-
    base(CR, V),
    RFC = label('Not-RFC-or-WIP', ok(_)).

base(CR, V) :-
    gerrit:max_with_block(-2, 2, 'Code-Review', CR),
    gerrit:max_with_block(-2, 2, 'Verified', V).
