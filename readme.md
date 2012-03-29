Abstract: A mathematical knot is simply a closed curve in three-space. Classifying open knots, or knots that have not been closed, is a relatively unexplored area of knot theory. In this note, we report on our study of open random walks of varying length, creating a collection of open knots. Following the strategy of Millett, Dobay and Stasiak, an open knot is closed by connecting its two open endpoints to a third point, lying on a large sphere that encloses the random walk deeply within its interior. The resulting polygonal knot can be analyzed and its knot type determined, up to the indetermincy of standard knot invariants, using the HOMFLY polynomial. With many closure points uniformly distributed on the large sphere, a statistical distribution of knot types is created for each open knot. We use this method to continue the exploration of the knottedness of linear random walks...

walker.f - Creates a random walk
knotter.center.f - Closes a set of random walks using the walk's centerpoint (based on radius of gyration), and analyses the resulting closed path using HOMFLY algorithm for classifying mathematical knots.
knotter.f - Closes the path using a random point on an enclosing sphere (I think)