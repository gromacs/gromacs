submit_rule(submit(CR, V)) :-
    base(CR, V),
    gerrit:commit_message_matches('(?!\\\[WIP\\\])|(?!\\\[RFC\\\])'),
    !.

submit_rule(submit(CR, V)) :-
    base(CR, V),
    label('Commit-Message-does-not-include-[WIP]-or-[RFC]', need(_)).

base(CR, V) :-
    gerrit:max_with_block(-2, 2, 'Code-Review', CR),
    gerrit:max_with_block(-2, 2, 'Verified', V).
