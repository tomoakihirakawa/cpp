#ifndef mooring_H
#define mooring_H

#include <array>

struct MooringLine {
  private:
   std::array<double, 3> anchorPosition;   // Anchor point in global coordinates
   std::array<double, 3> attachmentPoint;  // Attachment point on the floating body
   float lineLength;                       // Unstressed length of the mooring line
   float lineDiameter;                     // Diameter of the line
   Material material;                      // Material of the line
   float currentTension;                   // Current tension in the line

  public:
   // Constructor
   MooringLine(std::array<double, 3> anchor,
               std::array<double, 3> attach,
               float length,
               float diameter,
               Material mat);

   // Member functions
   void CalculateTension(const std::array<double, 3>& bodyPosition);                   // Calculate tension based on the floating body position
   std::array<double, 3> ApplyForceToBody(const std::array<double, 3>& bodyPosition);  // Apply force to floating body
   void Update(const std::array<double, 3>& bodyPosition);                             // Update line status
   std::array<double, 3> GetForceComponents();                                         // Get components of force
   bool CheckBreakingCondition();                                                      // Check if the line would break
   void CalculateWaveForce(const Wave& waveData);                                      // Calculate wave-induced forces, if applicable
   void LogStatus();                                                                   // Log current status of the line
};

#endif