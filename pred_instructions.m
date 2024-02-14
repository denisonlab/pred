function pred_instructions(window, p, devNum, data)


nInstruct = 7;
instructText{1} = sprintf('Slide 1 out of %d\n\nIn this experiment, you will make judgments about visual information.',nInstruct);
instructText{2} = sprintf('Slide 2 out of %d\n\nThe experiment has two parts. Today, you will complete one session of Part 1 and one session of Part 2. Part 1 will take approximately 10 min, and Part 2 will take approximately 45 min. We''ll take a look at instructions for Part 1, then show the instructions for Part 2 after you have completed Part 1.',nInstruct);
instructText{3} = sprintf('Slide 3 out of %d\n\nEach trial will begin with a central white dot.',nInstruct);
instructText{4} = sprintf('Slide 4 out of %d\n\nThen, a grating, tilted counterclockwise or clockwise from vertical, may or may not appear.',nInstruct);
instructText{5} = sprintf('Slide 5 out of %d\n\nWhen the dot turns from white to green, you will have to make a judgment about whether or not you saw a central grating. At the same time, you will also need to report whether the grating is counterclockwise or clockwise from vertical, regardless of if you feel you saw a grating or not.\n1→ see, counterclockwise\n2→ did not see, guess counterclockwise\n9→ did not see, guess clockwise\n0→ see, clockwise',nInstruct);
instructText{6} = sprintf('Slide 6 out of %d\n\nA good time to blink or look away from the dot is when the central dot turns green and you prepare your response. At all other times during the trial, it''s important to keep your gaze fixed on the dot so as to not miss a potential grating. In particular, resist the temptation to move your gaze toward the periphery. The eye tracker will record your eye movements throughout the trial.',nInstruct);
instructText{7} = sprintf('Slide 7 out of %d\n\nNext, we will look at a sample trial.',nInstruct);

waitText = sprintf('________\n\n Press space to continue');
endText = sprintf('________\n\n Press space to end');


for i = 1:nInstruct
    % Instruction text
    DrawFormattedText(window, instructText{i}, round(p.instruct.xPos*p.ppd), round(p.instruct.yPosUpper*p.ppd), [1 1 1]*s.white, p.instruct.wrapAt)
    timeInstruct = Screen('Flip', window);
    DrawFormattedText(window, instructText{i}, round(p.instruct.xPos*p.ppd), round(p.instruct.yPosUpper*p.ppd), [1 1 1]*s.white, p.instruct.wrapAt);
    % Continue text (after delay)
    DrawFormattedText(window, waitText, round(p.instruct.xPos*p.ppd), s.rect(4)-round(p.instruct.yPosLower*p.ppd), [1 1 1]*s.white);
    Screen('Flip', window, timeInstruct + p.instruct.waitTime);
    [secs, keyCode] = KbWait(devNum);
    responseKey = find(keyCode);
end

end
