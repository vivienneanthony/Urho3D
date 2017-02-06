//
// Copyright (c) 2008-2017 the Urho3D project.
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
#include "../Graphics/TerrainFace.h"
#include "../Graphics/Terrain.h"
#include "../Graphics/VertexBuffer.h"
#include "../IO/Log.h"
#include "../Scene/Node.h"

#include "../DebugNew.h"

namespace Urho3D
{

static const float LOD_CONSTANT = 1.0f / 150.0f;

extern const char* GEOMETRY_CATEGORY;


// Create each face detached
TerrainFace::TerrainFace(Context* context) :
    Drawable(context, DRAWABLE_GEOMETRY),
    geometry_(new Geometry(context)),
    maxLodGeometry_(new Geometry(context)),
    occlusionGeometry_(new Geometry(context)),
    vertexBuffer_(new VertexBuffer(context)),
    coordinates_(IntVector2::ZERO),
    m_LodLevel(0)
{
    // Create new geometry
    geometry_->SetVertexBuffer(0, vertexBuffer_);
    maxLodGeometry_->SetVertexBuffer(0, vertexBuffer_);
    occlusionGeometry_->SetVertexBuffer(0, vertexBuffer_);

    SetEnabled(true);
}


// Create each face detached
void TerrainFace::SetTerrainFace(unsigned int index, unsigned maxLod, Terrain * parentTerrain )
{
    m_Index=index;
    m_MaxLod=maxLod;
    m_ParentTerrain=parentTerrain;

    SetEnabled(true);
}


TerrainFace::~TerrainFace()
{
}

void TerrainFace::RegisterObject(Context* context)
{
    context->RegisterFactory<TerrainFace>();
}

void TerrainFace::ProcessRayQuery(const RayOctreeQuery& query, PODVector<RayQueryResult>& results)
{

}

void TerrainFace::UpdateBatches(const FrameInfo& frame)
{

}

void TerrainFace::UpdateGeometry(const FrameInfo& frame)
{

}

UpdateGeometryType TerrainFace::GetUpdateGeometryType()
{
    // Because there is a latency in starting worker thread updates, and the update of terrain patch LOD should not take
    // much time, always update in the main thread
    return UPDATE_MAIN_THREAD;
}

Geometry* TerrainFace::GetLodGeometry(unsigned batchIndex, unsigned level)
{
    if (!level)
        return maxLodGeometry_;
    else
        return geometry_;
}

unsigned TerrainFace::GetNumOccluderTriangles()
{
    // Check that the material is suitable for occlusion (default material always is)
    Material* mat = batches_[0].material_;
    if (mat && !mat->GetOcclusion())
        return 0;
    else
        return occlusionGeometry_->GetIndexCount() / 3;
}

bool TerrainFace::DrawOcclusion(OcclusionBuffer* buffer)
{
    const unsigned char* vertexData;
    unsigned vertexSize;
    const unsigned char* indexData;
    unsigned indexSize;
    const PODVector<VertexElement>* elements;

    geometry_->GetRawData(vertexData, vertexSize, indexData, indexSize, elements);
    // Check for valid geometry data
    if (!vertexData || !elements || VertexBuffer::GetElementOffset(*elements, TYPE_VECTOR3, SEM_POSITION) != 0)
    {
        return false;
    }

    // Draw and check for running out of triangles
    buffer->AddTriangles(node_->GetWorldTransform(), vertexData, vertexSize, geometry_->GetVertexStart(),
                         geometry_->GetVertexCount());


}

void TerrainFace::DrawDebugGeometry(DebugRenderer* debug, bool depthTest)
{

}

void TerrainFace::SetOwner(Terrain* terrain)
{

}

void TerrainFace::SetNeighbors(TerrainFace* north, TerrainFace* south, TerrainFace* west, TerrainFace* east)
{
    //north_ = north;
    //south_ = south;
    //west_ = west;
    //east_ = east;
}

void TerrainFace::SetMaterial(Material* material)
{

}

void TerrainFace::SetBoundingBox(const BoundingBox& box)
{

}

void TerrainFace::SetCoordinates(const IntVector2& coordinates)
{
    coordinates_ = coordinates;
}

void TerrainFace::ResetLod()
{
}

Geometry* TerrainFace::GetGeometry() const
{

}

Geometry* TerrainFace::GetMaxLodGeometry() const
{
    return maxLodGeometry_;
}

Geometry* TerrainFace::GetOcclusionGeometry() const
{
    return occlusionGeometry_;
}

VertexBuffer* TerrainFace::GetVertexBuffer() const
{
    return vertexBuffer_;
}

Terrain* TerrainFace::GetOwner() const
{
    return owner_;
}

void TerrainFace::OnWorldBoundingBoxUpdate()
{
}

unsigned TerrainFace::GetCorrectedLodLevel(unsigned lodLevel)
{
    return lodLevel;
}

void TerrainFace::BuildFace(Vector3 centre, Vector3 dx, Vector3 dy)
{
    // Define the grid size.  i.e. the number of vertices in a grid (16x16)
    unsigned int PatchSize = 16;
    unsigned int PatchSizePlus = 17; // PatchSize + 1

    // Set buffer
    vertexBuffer_->SetSize(PatchSize*PatchSize, MASK_POSITION | MASK_NORMAL | MASK_TEXCOORD1,true);

    Vector<float> _vertexData;

    // Didn't work
    float* tempVertexBuffer = (float*)vertexBuffer_->Lock(0, vertexBuffer_->GetVertexCount());

    int xPos=0;
    int yPos=0;

    // loop through and create a grid of vertices.
    for (int u = 0; u <= PatchSize; u++)
    {
        for (int v = 0; v <= PatchSize; v++)
        {

            // Create the vertex grid around the centre of thecube face (which is passed into the function as Vector3 centre).
            Vector3 tempPosition = centre + (dx / PatchSize) * (v - PatchSize / 2) + (dy / PatchSize) * (u - PatchSize / 2);

            // This is where we would define the height of the vertex.
            float fheight = 0;

            // Project the vertex onto the sphere, taking into consideration the height of the
            // vertex and the radius of the planet.  By specifying 0 as the height, we will
            // get a 'perfectly' round planet/sphere.
            Vector3 newPosition = SurfaceVectorToCoordinates(tempPosition, m_ParentTerrain->GetWorldRadius(), fheight);

            // copy vertex
            _vertexData.Push(newPosition.x_);
            _vertexData.Push(newPosition.y_);
            _vertexData.Push(newPosition.z_);

            // Normal
            _vertexData.Push(1.0f);
            _vertexData.Push(1.0f);
            _vertexData.Push(1.0f); //should point to parent center)

            // Normal
            Vector2 texCoord((float)u/ (float)( PatchSize - 1), 1.0f - (float)v / (float)( PatchSize  - 1));
            _vertexData.Push(texCoord.x_);
            _vertexData.Push(texCoord.y_);
        }
    }

    vertexBuffer_->SetData((float *) tempVertexBuffer);

    vertexBuffer_->Unlock();
}

/// <summary>
/// Transforms a vector from the surface of a cube, onto the surface of a sphere.
/// </summary>
/// surfacePos is the vertex position on the cube.
/// radius is the radius of the planet.
/// height is the height of the terrain at this position.
Vector3 TerrainFace::SurfaceVectorToCoordinates(Vector3 surfacePos, float radius, float height)
{

// Create a return veriable.
    Vector3 loReturnData = surfacePos;

// Get a unit vector ( this will 'point' in the correct direction, from (0,0,0) to
// the position of the vertex on the sphere ).
    loReturnData.Normalize();

// Add the planet radius and the height of the vertex, and return the vector.
    return loReturnData *(radius + height);
}


}
