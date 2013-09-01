submit_rule(submit(CR, V, NoSubmitTags)) :-
    base(CR, V),
    gerrit:commit_message_matches('(?!\\\[WIP\\\])|(?!\\\[RFC\\\])'),
    NoSubmitTags = label(_, ok(_)),
    !.

submit_rule(submit(CR, V, NoSubmitTags)) :-
    base(CR, V),
    NoSubmitTags = label('Commit-Message-does-not-include-[WIP]-or-[RFC]', need(_)).

base(CR, V) :-
    gerrit:max_with_block(-2, 2, 'Code-Review', CR),
    gerrit:max_with_block(-2, 2, 'Verified', V).
