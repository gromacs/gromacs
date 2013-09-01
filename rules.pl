submit_rule(submit(WIP)) :-
    WIP = label('Commit-Message-does-not-include-[WIP]', need(_)).

submit_rule(submit(RFC)) :-
    RFC = label('Commit-Message-does-not-include-[RFC]', need(_)).

submit_rule(submit(CR, V)) :-
    gerrit:commit_message_matches('(?!\\\[WIP\\\])'),
    gerrit:commit_message_matches('(?!\\\[RFC\\\])'),
    gerrit:max_with_block(-2, 2, 'Code-Review', CR),
    gerrit:max_with_block(-2, 2, 'Verified', V).
