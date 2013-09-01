submit_rule(submit(V, CR)) :-
    gerrit:max_with_block(-2, 2, 'Code-Review', CR),
    gerrit:max_with_block(-2, 2, 'Verified', V).

submit_rule(submit()) :-
    not(gerrit:commit_message_matches('\\\[WIP\\\]')),
    WIP = label('Commit-Message-does-not-include-[WIP]', need(_)).

submit_rule(submit()) :-
    not(gerrit:commit_message_matches('\\\[RFC\\\]')),
    RFC = label('Commit-Message-does-not-include-[RFC]', need(_)).
