//
// Copyright (c) 2008-2016 the Urho3D project.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//

#include "../Precompiled.h"

#include "../Core/Context.h"
#include "../Graphics/Camera.h"
#include "../Graphics/DebugRenderer.h"
#include "../Graphics/Geometry.h"
#include "../Graphics/IndexBuffer.h"
#include "../Graphics/Material.h"
#include "../Graphics/OcclusionBuffer.h"
#include "../Graphics/OctreeQuery.h"
#include "../Graphics/Terrain.h"
#include "../Graphics/TerrainPatch.h"
#include "../Graphics/VertexBuffer.h"
#include "../IO/Log.h"
#include "../Scene/Node.h"

#include "../DebugNew.h"

namespace Urho3D
{

static const float LOD_CONSTANT = 1.0f / 150.0f;

extern const char* GEOMETRY_CATEGORY;

TerrainPatch::TerrainPatch(Context* context) :
    Drawable(context, DRAWABLE_GEOMETRY),
    geometry_(new Geometry(context)),
    maxLodGeometry_(new Geometry(context)),
    occlusionGeometry_(new Geometry(context)),
    vertexBuffer_(new VertexBuffer(context)),
    coordinates_(IntVector2::ZERO),
    lodLevel_(0)
{
    geometry_->SetVertexBuffer(0, vertexBuffer_);
    maxLodGeometry_->SetVertexBuffer(0, vertexBuffer_);
    occlusionGeometry_->SetVertexBuffer(0, vertexBuffer_);

    batches_.Resize(1);
    batches_[0].geometry_ = geometry_;
    batches_[0].geometryType_ = GEOM_STATIC_NOINSTANCING;
}


// New Construct
TerrainPatch::TerrainPatch(Context* context,unsigned int index, unsigned int maxLOD, Terrain * parentTerrain):
    Drawable(GetContext(), DRAWABLE_GEOMETRY),
    geometry_(new Geometry(context)),
    maxLodGeometry_(new Geometry(context)),
    occlusionGeometry_(new Geometry(context)),
    vertexBuffer_(new VertexBuffer(context)),
    coordinates_(IntVector2::ZERO),
    lodLevel_(0),
    m_MaxLod(maxLOD),
    m_Index(index),
    m_ParentTerrain(parentTerrain),
    indexBuffer_(new IndexBuffer(context))
{
    geometry_->SetVertexBuffer(0, vertexBuffer_);
    maxLodGeometry_->SetVertexBuffer(0, vertexBuffer_);
    occlusionGeometry_->SetVertexBuffer(0, vertexBuffer_);

    batches_.Resize(1);
    batches_[0].geometry_ = geometry_;
    batches_[0].geometryType_ = GEOM_STATIC_NOINSTANCING;

}


TerrainPatch::~TerrainPatch()
{
}

void TerrainPatch::RegisterObject(Context* context)
{
    context->RegisterFactory<TerrainPatch>();
}

void TerrainPatch::ProcessRayQuery(const RayOctreeQuery& query, PODVector<RayQueryResult>& results)
{
    /*  RayQueryLevel level = query.level_;

      switch (level)
      {
      case RAY_AABB:
          Drawable::ProcessRayQuery(query, results);
          break;

      case RAY_OBB:
      case RAY_TRIANGLE:
          {
              Matrix3x4 inverse(node_->GetWorldTransform().Inverse());
              Ray localRay = query.ray_.Transformed(inverse);
              float distance = localRay.HitDistance(boundingBox_);
              Vector3 normal = -query.ray_.direction_;

              if (level == RAY_TRIANGLE && distance < query.maxDistance_)
              {
                  Vector3 geometryNormal;
                  distance = geometry_->GetHitDistance(localRay, &geometryNormal);
                  normal = (node_->GetWorldTransform() * Vector4(geometryNormal, 0.0f)).Normalized();
              }

              if (distance < query.maxDistance_)
              {
                  RayQueryResult result;
                  result.position_ = query.ray_.origin_ + distance * query.ray_.direction_;
                  result.normal_ = normal;
                  result.distance_ = distance;
                  result.drawable_ = this;
                  result.node_ = node_;
                  result.subObject_ = M_MAX_UNSIGNED;
                  results.Push(result);
              }
          }
          break;

      case RAY_TRIANGLE_UV:
          URHO3D_LOGWARNING("RAY_TRIANGLE_UV query level is not supported for TerrainPatch component");
          break;
      }*/
}

void TerrainPatch::UpdateBatches(const FrameInfo& frame)
{
    /*const Matrix3x4& worldTransform = node_->GetWorldTransform();
    distance_ = frame.camera_->GetDistance(GetWorldBoundingBox().Center());

    float scale = worldTransform.Scale().DotProduct(DOT_SCALE);
    lodDistance_ = frame.camera_->GetLodDistance(distance_, scale, lodBias_);

    batches_[0].distance_ = distance_;
    batches_[0].worldTransform_ = &worldTransform;

    unsigned newLodLevel = 0;
    for (unsigned i = 0; i < lodErrors_.Size(); ++i)
    {
        if (lodErrors_[i] / lodDistance_ > LOD_CONSTANT)
            break;
        else
            newLodLevel = i;
    }

    lodLevel_ = GetCorrectedLodLevel(newLodLevel);*/
}

void TerrainPatch::UpdateGeometry(const FrameInfo& frame)
{
    /* if (vertexBuffer_->IsDataLost())
     {
         if (owner_)
             owner_->CreatePatchGeometry(this);
         else
             vertexBuffer_->ClearDataLost();
     }

     if (owner_)
         owner_->UpdatePatchLod(this);*/
}

UpdateGeometryType TerrainPatch::GetUpdateGeometryType()
{
    // Because there is a latency in starting worker thread updates, and the update of terrain patch LOD should not take
    // much time, always update in the main thread
    return UPDATE_MAIN_THREAD;
}

Geometry* TerrainPatch::GetLodGeometry(unsigned batchIndex, unsigned level)
{
    if (!level)
        return maxLodGeometry_;
    else
        return geometry_;
}

unsigned TerrainPatch::GetNumOccluderTriangles()
{
    // Check that the material is suitable for occlusion (default material always is)
    Material* mat = batches_[0].material_;
    if (mat && !mat->GetOcclusion())
        return 0;
    else
        return occlusionGeometry_->GetIndexCount();
}

bool TerrainPatch::DrawOcclusion(OcclusionBuffer* buffer)
{
    // Check that the material is suitable for occlusion (default material always is) and set culling mode
    Material* material = batches_[0].material_;
    if (material)
    {
        if (!material->GetOcclusion())
            return true;
        buffer->SetCullMode(material->GetCullMode());
    }
    else
        buffer->SetCullMode(CULL_CCW);

    const unsigned char* vertexData;
    unsigned vertexSize;
    const unsigned char* indexData;
    unsigned indexSize;
    const PODVector<VertexElement>* elements;

    occlusionGeometry_->GetRawData(vertexData, vertexSize, indexData, indexSize, elements);
    // Check for valid geometry data
    if (!vertexData || !indexData || !elements || VertexBuffer::GetElementOffset(*elements, TYPE_VECTOR3, SEM_POSITION) != 0)
        return false;

    // DDraw and check for running out of triangles
    return buffer->AddTriangles(node_->GetWorldTransform(), vertexData, vertexSize, indexData, indexSize, occlusionGeometry_->GetIndexStart(),
                                occlusionGeometry_->GetIndexCount());
}

void TerrainPatch::DrawDebugGeometry(DebugRenderer* debug, bool depthTest)
{
    // Intentionally no operation
}

void TerrainPatch::SetOwner(Terrain* terrain)
{
    owner_ = terrain;
}

void TerrainPatch::SetNeighbors(TerrainPatch* north, TerrainPatch* south, TerrainPatch* west, TerrainPatch* east)
{
    north_ = north;
    south_ = south;
    west_ = west;
    east_ = east;
}

void TerrainPatch::SetMaterial(Material* material)
{
    batches_[0].material_ = material;
}

void TerrainPatch::SetBoundingBox(const BoundingBox& box)
{
    boundingBox_ = box;
    OnMarkedDirty(node_);
}

void TerrainPatch::SetCoordinates(const IntVector2& coordinates)
{
    coordinates_ = coordinates;
}

void TerrainPatch::ResetLod()
{
    lodLevel_ = 0;
}

Geometry* TerrainPatch::GetGeometry() const
{
    return geometry_;
}

Geometry* TerrainPatch::GetMaxLodGeometry() const
{
    return maxLodGeometry_;
}

Geometry* TerrainPatch::GetOcclusionGeometry() const
{
    return occlusionGeometry_;
}

VertexBuffer* TerrainPatch::GetVertexBuffer() const
{
    return vertexBuffer_;
}

Terrain* TerrainPatch::GetOwner() const
{
    return owner_;
}

void TerrainPatch::OnWorldBoundingBoxUpdate()
{
    worldBoundingBox_ = boundingBox_.Transformed(node_->GetWorldTransform());
}

unsigned TerrainPatch::GetCorrectedLodLevel(unsigned lodLevel)
{
    if (north_)
        lodLevel = Min(lodLevel, north_->GetLodLevel() + 1);
    if (south_)
        lodLevel = Min(lodLevel, south_->GetLodLevel() + 1);
    if (west_)
        lodLevel = Min(lodLevel, west_->GetLodLevel() + 1);
    if (east_)
        lodLevel = Min(lodLevel, east_->GetLodLevel() + 1);

    return lodLevel;
}



/// <summary
/// Builds the parent node information for this quadtree face.  The actual terrain mesh is generated in this function.
/// Note that for child nodes, there is a separate function, BuildChildNode().
/// </summary>
/// game is the current instance of the XNA game class.
/// centre is the centre of the planet face (on one of the cube edges).
/// dx is the direction of the X-Axis in local 'terrain space'.
/// dy is the direction of the Y-Axis in local 'terrain space'.
void TerrainPatch::BuildParentNode(Vector3 centre, Vector3 dx, Vector3 dy)
{
    // Define the grid size.  i.e. the number of vertices in a grid (16x16)
    unsigned int GridSize = 16;
    unsigned int GridSizePlus = 17; // GridSize + 1

    // Create a new buffer for data
    unsigned int * rawDataVertex =  (unsigned int *)new float[GridSize * GridSize * sizeof(Vector3)];

    // Set buffer to this
    vertexBuffer_->SetData(rawDataVertex);

    // Set buffer size - Set the normal and size
    vertexBuffer_->SetSize(GridSize*GridSize, MASK_POSITION | MASK_NORMAL | MASK_TEXCOORD1,true);

    float* newDataSize = (float*)vertexBuffer_->Lock(0, vertexBuffer_->GetVertexCount());


    BoundingBox box;

    // loop through and create a grid of vertices.
    for (int u = 0; u <= GridSize; u++)
    {
        for (int v = 0; v <= GridSize; v++)
        {
            // Create the vertex grid around the centre of thecube face (which is passed into the function as Vector3 centre).
            Vector3 tempPosition = centre + (dx / GridSize) * (v - GridSize / 2) + (dy / GridSize) * (u - GridSize / 2);

            // This is where we would define the height of the vertex.
            float lfHeight = 0.0f;

            // Project the vertex onto the sphere, taking into consideration the height of the
            // vertex and the radius of the planet.  By specifying 0 as the height, we will
            // get a 'perfectly' round planet/sphere.
            unsigned int position = GridSizePlus * u + v;

            Vector3 pos = SurfaceVectorToCoordinates(tempPosition, 1.0f, lfHeight);

            // set position
            newDataSize[position] = pos.x_;
            newDataSize[position+1] = pos.y_;
            newDataSize[position+2] = pos.z_;

            box.Merge(pos);

        }
    }



    geometry_->SetIndexBuffer(indexBuffer_);
    //geometry_->SetRawVertexData(vertexBuffer_, MASK_POSITION);
    maxLodGeometry_->SetIndexBuffer(indexBuffer_);
    //maxLodGeometry_->SetRawVertexData(vertexBuffer_, MASK_POSITION);
    occlusionGeometry_->SetIndexBuffer(indexBuffer_);
   //occlusionGeometry_->SetRawVertexData(vertexBuffer_, MASK_POSITION);


    // Create the shared index data
    CreateIndexData();

    vertexBuffer_->Unlock();

    //CreateIndexData();

}

/// <summary>
/// Transforms a vector from the surface of a cube, onto the surface of a sphere.
/// </summary>
/// surfacePos is the vertex position on the cube.
/// radius is the radius of the planet.
/// height is the height of the terrain at this position.
Vector3 TerrainPatch::SurfaceVectorToCoordinates(Vector3 surfacePos, float radius, float height)
{

    // Create a return veriable.
    Vector3 loReturnData = surfacePos;

    // Get a unit vector ( this will 'point' in the correct direction, from (0,0,0) to
    // the position of the vertex on the sphere ).
    loReturnData.Normalize();

    // Add the planet radius and the height of the vertex, and return the vector.
    return loReturnData *(radius + height);
}


void TerrainPatch::CreateIndexData()
{
    PODVector<unsigned short> indices;
    drawRanges_.Clear();

    /* Build index data for each LOD level. Each LOD level except the lowest can stitch to the next lower LOD from the edges:
       north, south, west, east, or any combination of them, requiring 16 different versions of each LOD level's index data
       Normal edge:     Stitched edge:
       +----+----+      +---------+
       |\   |\   |      |\       /|
       | \  | \  |      | \     / |
       |  \ |  \ |      |  \   /  |
       |   \|   \|      |   \ /   |
       +----+----+      +----+----+
    */
    unsigned indexStart = indices.Size();

    int zStart = 0;
    int xStart = 0;
    int zEnd = 16;
    int xEnd = 16;

    int skip = 1;

    unsigned row = (unsigned)17;


    // Build the main grid
    for (int z = zStart; z < zEnd; z += skip)
    {
        for (int x = xStart; x < xEnd; x += skip)
        {
            indices.Push((unsigned short)((z + skip) * row + x));
            indices.Push((unsigned short)(z * row + x + skip));
            indices.Push((unsigned short)(z * row + x));
            indices.Push((unsigned short)((z + skip) * row + x));
            indices.Push((unsigned short)((z + skip) * row + x + skip));
            indices.Push((unsigned short)(z * row + x + skip));
        }
    }

    drawRanges_.Push(MakePair(indexStart, indices.Size() - indexStart));

    indexBuffer_->SetSize(indices.Size(), false);
    indexBuffer_->SetData(&indices[0]);

    geometry_->SetDrawRange(TRIANGLE_LIST, indexStart, indices.Size() - indexStart, false);
    maxLodGeometry_->SetDrawRange(TRIANGLE_LIST, indexStart, indices.Size() - indexStart, false);
    occlusionGeometry_->SetDrawRange(TRIANGLE_LIST, indexStart, indices.Size() - indexStart, false);


}



}




