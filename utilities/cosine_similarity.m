function similarity = cosine_similarity(v1, v2)
%% Function to get cosine similarity of two vectors
% Should be bounded on [-1 1]

v1v2 = dot(v1, v2);
v1_norm = norm(v1);
v2_norm = norm(v2);

v1v2_norm = v1_norm*v2_norm;

similarity = v1v2/v1v2_norm;
end