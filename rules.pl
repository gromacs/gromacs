submit_rule(submit(V, CR)) :-
    gerrit:max_with_block(-2, 2, 'Code-Review', CR),
    gerrit:max_with_block(-2, 2, 'Verified', V).

submit_rule(submit(WIP)) :-
    gerrit:commit_message_matches('WIP'),
    WIP = label('Commit-Message-includes-WIP', reject(_)).
