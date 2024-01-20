% Version 1: 1.18.2024
% A simple script that defines a class used to create CelestialBodies to be used in other scripts. Right now, aside from the bare necesities (Velocity and Position calculations, other parts of the script are non-functioning)

classdef CelestialBody < handle
    properties 
        initialPosition {mustBeNumeric};
        initialVelocity {mustBeNumeric};
        Position        {mustBeNumeric};
        Velocity        {mustBeNumeric};
        Mass            {mustBeNumeric};
        Radius          {mustBeNumeric};
        ID;
    end
    methods
        function CB = CelestialBody(initposition, initvelocity, Mass, Radius) % Constructor
            CB.initialPosition = initposition;
            CB.initialVelocity = initvelocity;
            CB.Mass = Mass;
            CB.Radius = Radius;
            CB.Position = CB.initialPosition;
            CB.Velocity = CB.initialVelocity;
            Index(CB);
        end

        % Indexing is separate right now since it doesn't work. The idea was to use something like java.rmi.server.UID(),
        % but ordering as an integer starting at 1 to go along when all celestial bodies are put together,
        % But managing that (without creating another class or subfunction) seems to be difficult
        % Will return to this later

        function CelestialBodies = Index(CB) 
            CB.ID = randi(400); 
        end

        % Update Velocity, following common computation laws. One great thing about MATLAB is 
        % the fact that every variable called is immediately a double and not a float, as that
        % would cause problems. Of course, C# and C++ allow you to declare variable types, 
        % but that also leads to some to think that a float would be enough for every situation. 
        function UpdateVelocity(CB, CelestialBodies, dt)
            gravConst = 6.67430e-20;
            acceleration = zeros(1,3);
            for i = 1:length(CelestialBodies)
                if ~isequal(CB, CelestialBodies(i))
                    distance = CelestialBodies(i).Position - CB.Position;
                    acceleration = acceleration + gravConst * CelestialBodies(i).Mass * distance / norm(distance)^3;
                end
            end
            CB.Velocity = CB.Velocity + acceleration * dt;
        end

        % FTWDK: Updating velocity and position together is a big no-no.
        % In real life, Velocity and Position are continuously changing together. For
        % simulation, every body has to know its velocity changing before moving, otherwise
        % you will get a leapfrog-like effect, which will not be true to reality.
        % (Many people make this mistake early on, but not me, I'm Him)
        function UpdatePosition(CB,dt)
            CB.Position = CB.Position + CB.Velocity * dt;
        end
    end
end
