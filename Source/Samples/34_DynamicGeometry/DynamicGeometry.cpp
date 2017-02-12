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

#include <Urho3D/Core/CoreEvents.h>
#include <Urho3D/Core/Profiler.h>
#include <Urho3D/Engine/Engine.h>
#include <Urho3D/Graphics/Camera.h>
#include <Urho3D/Graphics/Geometry.h>
#include <Urho3D/Graphics/Graphics.h>
#include <Urho3D/Graphics/IndexBuffer.h>
#include <Urho3D/Graphics/Light.h>
#include <Urho3D/Graphics/Model.h>
#include <Urho3D/Graphics/Octree.h>
#include <Urho3D/Graphics/Renderer.h>
#include <Urho3D/Graphics/StaticModel.h>
#include <Urho3D/Graphics/VertexBuffer.h>
#include <Urho3D/Graphics/Zone.h>
#include <Urho3D/Input/Input.h>
#include <Urho3D/IO/Log.h>
#include <Urho3D/Resource/ResourceCache.h>
#include <Urho3D/Scene/Scene.h>
#include <Urho3D/UI/Font.h>
#include <Urho3D/UI/Text.h>
#include <Urho3D/UI/UI.h>

#include "DynamicGeometry.h"

#include <Urho3D/DebugNew.h>

URHO3D_DEFINE_APPLICATION_MAIN(DynamicGeometry)

DynamicGeometry::DynamicGeometry(Context* context) :
    Sample(context),
    animate_(true),
    time_(0.0f)
{
}

void DynamicGeometry::Start()
{
    // Execute base class startup
    Sample::Start();

    // Create the scene content
    CreateScene();

    // Create the UI content
    CreateInstructions();

    // Setup the viewport for displaying the scene
    SetupViewport();

    // Hook up to the frame update events
    SubscribeToEvents();

    // Set the mouse mode to use in the sample
    Sample::InitMouseMode(MM_RELATIVE);
}

void DynamicGeometry::CreateScene()
{
    ResourceCache* cache = GetSubsystem<ResourceCache>();

    scene_ = new Scene(context_);

    // Create the Octree component to the scene so that drawable objects can be rendered. Use default volume
    // (-1000, -1000, -1000) to (1000, 1000, 1000)
    scene_->CreateComponent<Octree>();

    // Create a Zone for ambient light & fog control
    Node* zoneNode = scene_->CreateChild("Zone");
    Zone* zone = zoneNode->CreateComponent<Zone>();
    zone->SetBoundingBox(BoundingBox(-1000.0f, 1000.0f));
    zone->SetFogColor(Color(0.2f, 0.2f, 0.2f));
    zone->SetFogStart(200.0f);
    zone->SetFogEnd(300.0f);

    // Create a directional light
    Node* lightNode = scene_->CreateChild("DirectionalLight");
    lightNode->SetDirection(Vector3(-0.6f, -1.0f, -0.8f)); // The direction vector does not need to be normalized
    Light* light = lightNode->CreateComponent<Light>();
    light->SetLightType(LIGHT_DIRECTIONAL);
    light->SetColor(Color(0.4f, 1.0f, 0.4f));
    light->SetSpecularIntensity(1.5f);

    // Get the original model and its unmodified vertices, which are used as source data for the animation
    Model* originalModel = cache->GetResource<Model>("Models/Box.mdl");
    if (!originalModel)
    {
        URHO3D_LOGERROR("Model not found, cannot initialize example scene");
        return;
    }



    // Finally create one model (pyramid shape) and a StaticModel to display it from scratch
    // Note: there are duplicated vertices to enable face normals. We will calculate normals programmatically
    unsigned numVertices = 18;


    // set chunks or patches
    int PatchSize=8;


    // Set Cube size
    float CubeSize=10;;

    // Get Cube Center - Create Defaults
    Vector3 center = Vector3::ZERO;
    Vector3 direction_x = Vector3::ZERO;
    Vector3 direction_y = Vector3::ZERO;


    // Create buffer
    Vector<Vector3> vertexData;
    Vector<unsigned short> indexData;

    // Create a index
    unsigned int index = 0;

    for(unsigned int face=0; face< 6; face++)
    {

        // Get Cube Center
        center = ((float)CubeSize/2)*TerrainFaceCoordinate[face];

        // Calculate direction based of x
        switch(face)
        {
        case Pos_X:
            direction_x = CubeSize*Vector3(0,1,0);
            direction_y = CubeSize*Vector3(0,0,1);
            break;
        case Neg_X:
            direction_x = CubeSize*Vector3(0,1,0);
            direction_y = CubeSize*Vector3(0,0,1);
            break;
        case Pos_Y:
            direction_x = CubeSize*Vector3(1,0,0);
            direction_y = CubeSize*Vector3(0,0,1);
            break;
        case Neg_Y:
            direction_x = CubeSize*Vector3(1,0,0);
            direction_y = CubeSize*Vector3(0,0,1);
            break;
        case Pos_Z:
            direction_x = CubeSize*Vector3(1,0,0);
            direction_y = CubeSize*Vector3(0,1,0);
            break;
        case Neg_Z:
            direction_x = CubeSize*Vector3(1,0,0);
            direction_y = CubeSize*Vector3(0,1,0);
            break;
        }

        // loop through and create a grid of vertices. // do not draw edge
        for (int u = 0; u < PatchSize; u++)
        {
            for (int v = 0; v < PatchSize; v++)
            {
                // Calculate patch size
                Vector3 x0=(direction_x / PatchSize) * (v- PatchSize / 2);
                Vector3 y0=(direction_y / PatchSize) * (u- PatchSize / 2);

                // Create the vertex grid around the center of thecube face (which is passed into the function as Vector3 center).
                Vector3 v0= center+ x0+y0;

                // Calculate patch size
                Vector3 x1=(direction_x / PatchSize) * (v - PatchSize / 2);
                Vector3 y1=(direction_y / PatchSize) * (u+1- PatchSize / 2);

                // Create the vertex grid around the center of thecube face (which is passed into the function as Vector3 center).
                Vector3 v1= center+ x1+y1;

                // Calculate patch size
                Vector3 x2=(direction_x / PatchSize) * (v+1- PatchSize / 2);
                Vector3 y2=(direction_y / PatchSize) * (u- PatchSize / 2);

                // Create the vertex grid around the center of thecube face (which is passed into the function as Vector3 center).
                Vector3 v2= center+ x2+y2;

                // Calculate patch size
                Vector3 x3=(direction_x / PatchSize) * ((v+1)- PatchSize / 2);
                Vector3 y3=(direction_y / PatchSize) * ((u+1)- PatchSize / 2);

                // Create the vertex grid around the center of thecube face (which is passed into the function as Vector3 center).
                Vector3 v3=center+ x3+y3;

                Vector3 edge1 = v2-v1;
                Vector3 edge2 = v3-v1;

                Vector3 n1 = edge1.CrossProduct(edge2).Normalized();

                // order is wrong
                vertexData.Push(v0);
                vertexData.Push(n1);
                indexData.Push(index++);

                vertexData.Push(v1);
                vertexData.Push(n1);
                indexData.Push(index++);



                vertexData.Push(v2);
                vertexData.Push(n1);
                indexData.Push(index++);


                vertexData.Push(v1);
                vertexData.Push(n1);
                indexData.Push(index++);



                vertexData.Push(v2);
                vertexData.Push(n1);
                indexData.Push(index++);


                vertexData.Push(v3);
                vertexData.Push(n1);
                indexData.Push(index++);


            }
        }
    }
    numVertices = indexData.Size();

    SharedPtr<Model> fromScratchModel(new Model(context_));
    SharedPtr<VertexBuffer> vb(new VertexBuffer(context_));
    SharedPtr<IndexBuffer> ib(new IndexBuffer(context_));
    SharedPtr<Geometry> geom(new Geometry(context_));

    // Shadowed buffer needed for raycasts to work, and so that data can be automatically restored on device loss
    vb->SetShadowed(true);

    // We could use the "legacy" element bitmask to define elements for more compact code, but let's demonstrate
    // defining the vertex elements explicitly to allow any element types and order
    PODVector<VertexElement> elements;
    elements.Push(VertexElement(TYPE_VECTOR3, SEM_POSITION));
    elements.Push(VertexElement(TYPE_VECTOR3, SEM_NORMAL));


    vb->SetSize(numVertices, elements);
    vb->SetData(&vertexData[0]);

    ib->SetShadowed(true);
    ib->SetSize(numVertices, false);
    ib->SetData(&indexData[0]);

    geom->SetVertexBuffer(0, vb);
    geom->SetIndexBuffer(ib);
    geom->SetDrawRange(TRIANGLE_LIST, 0, numVertices);

    fromScratchModel->SetNumGeometries(1);
    fromScratchModel->SetGeometry(0, 0, geom);
    fromScratchModel->SetBoundingBox(BoundingBox(Vector3(-0.5f, -0.5f, -0.5f), Vector3(0.5f, 0.5f, 0.5f)));

    // Though not necessary to render, the vertex & index buffers must be listed in the model so that it can be saved properly
    Vector<SharedPtr<VertexBuffer> > vertexBuffers;
    Vector<SharedPtr<IndexBuffer> > indexBuffers;
    vertexBuffers.Push(vb);
    indexBuffers.Push(ib);

    // Morph ranges could also be not defined. Here we simply define a zero range (no morphing) for the vertex buffer
    PODVector<unsigned> morphRangeStarts;
    PODVector<unsigned> morphRangeCounts;
    morphRangeStarts.Push(0);
    morphRangeCounts.Push(0);
    fromScratchModel->SetVertexBuffers(vertexBuffers, morphRangeStarts, morphRangeCounts);
    fromScratchModel->SetIndexBuffers(indexBuffers);

    Node* node = scene_->CreateChild("FromScratchObject");
    node->SetPosition(Vector3(0.0f, 3.0f, 0.0f));
    StaticModel* object = node->CreateComponent<StaticModel>();
    object->SetModel(fromScratchModel);

    // Create the camera
    cameraNode_ = new Node(context_);
    cameraNode_->SetPosition(Vector3(0.0f, 2.0f, -20.0f));
    Camera* camera = cameraNode_->CreateComponent<Camera>();
    camera->SetFarClip(300.0f);
}

void DynamicGeometry::CreateInstructions()
{
    ResourceCache* cache = GetSubsystem<ResourceCache>();
    UI* ui = GetSubsystem<UI>();

    // Construct new Text object, set string to display and font to use
    Text* instructionText = ui->GetRoot()->CreateChild<Text>();
    instructionText->SetText(
        "Use WASD keys and mouse/touch to move\n"
        "Space to toggle animation"
    );
    instructionText->SetFont(cache->GetResource<Font>("Fonts/Anonymous Pro.ttf"), 15);
    // The text has multiple rows. Center them in relation to each other
    instructionText->SetTextAlignment(HA_CENTER);

    // Position the text relative to the screen center
    instructionText->SetHorizontalAlignment(HA_CENTER);
    instructionText->SetVerticalAlignment(VA_CENTER);
    instructionText->SetPosition(0, ui->GetRoot()->GetHeight() / 4);
}

void DynamicGeometry::SetupViewport()
{
    Renderer* renderer = GetSubsystem<Renderer>();

    // Set up a viewport to the Renderer subsystem so that the 3D scene can be seen
    SharedPtr<Viewport> viewport(new Viewport(context_, scene_, cameraNode_->GetComponent<Camera>()));
    renderer->SetViewport(0, viewport);
}

void DynamicGeometry::SubscribeToEvents()
{
    // Subscribe HandleUpdate() function for processing update events
    SubscribeToEvent(E_UPDATE, URHO3D_HANDLER(DynamicGeometry, HandleUpdate));
}

void DynamicGeometry::MoveCamera(float timeStep)
{
    // Do not move if the UI has a focused element (the console)
    if (GetSubsystem<UI>()->GetFocusElement())
        return;

    Input* input = GetSubsystem<Input>();

    // Movement speed as world units per second
    const float MOVE_SPEED = 20.0f;
    // Mouse sensitivity as degrees per pixel
    const float MOUSE_SENSITIVITY = 0.1f;

    // Use this frame's mouse motion to adjust camera node yaw and pitch. Clamp the pitch between -90 and 90 degrees
    IntVector2 mouseMove = input->GetMouseMove();
    yaw_ += MOUSE_SENSITIVITY * mouseMove.x_;
    pitch_ += MOUSE_SENSITIVITY * mouseMove.y_;
    pitch_ = Clamp(pitch_, -90.0f, 90.0f);

    // Construct new orientation for the camera scene node from yaw and pitch. Roll is fixed to zero
    cameraNode_->SetRotation(Quaternion(pitch_, yaw_, 0.0f));

    // Read WASD keys and move the camera scene node to the corresponding direction if they are pressed
    if (input->GetKeyDown(KEY_W))
        cameraNode_->Translate(Vector3::FORWARD * MOVE_SPEED * timeStep);
    if (input->GetKeyDown(KEY_S))
        cameraNode_->Translate(Vector3::BACK * MOVE_SPEED * timeStep);
    if (input->GetKeyDown(KEY_A))
        cameraNode_->Translate(Vector3::LEFT * MOVE_SPEED * timeStep);
    if (input->GetKeyDown(KEY_D))
        cameraNode_->Translate(Vector3::RIGHT * MOVE_SPEED * timeStep);
}

void DynamicGeometry::AnimateObjects(float timeStep)
{
    URHO3D_PROFILE(AnimateObjects);
    /*
        time_ += timeStep * 100.0f;

        // Repeat for each of the cloned vertex buffers
        for (unsigned i = 0; i < animatingBuffers_.Size(); ++i)
        {
            float startPhase = time_ + i * 30.0f;
            VertexBuffer* buffer = animatingBuffers_[i];

            // Lock the vertex buffer for update and rewrite positions with sine wave modulated ones
            // Cannot use discard lock as there is other data (normals, UVs) that we are not overwriting
            unsigned char* vertexData = (unsigned char*)buffer->Lock(0, buffer->GetVertexCount());
            if (vertexData)
            {
                unsigned vertexSize = buffer->GetVertexSize();
                unsigned numVertices = buffer->GetVertexCount();
                for (unsigned j = 0; j < numVertices; ++j)
                {
                    // If there are duplicate vertices, animate them in phase of the original
                    float phase = startPhase + vertexDuplicates_[j] * 10.0f;
                    Vector3& src = originalVertices_[j];
                    Vector3& dest = *reinterpret_cast<Vector3*>(vertexData + j * vertexSize);
                    dest.x_ = src.x_ * (1.0f + 0.1f * Sin(phase));
                    dest.y_ = src.y_ * (1.0f + 0.1f * Sin(phase + 60.0f));
                    dest.z_ = src.z_ * (1.0f + 0.1f * Sin(phase + 120.0f));
                }

                buffer->Unlock();
            }
        }
        */
}

void DynamicGeometry::HandleUpdate(StringHash eventType, VariantMap& eventData)
{
    using namespace Update;

    // Take the frame time step, which is stored as a float
    float timeStep = eventData[P_TIMESTEP].GetFloat();

    // Toggle animation with space
    Input* input = GetSubsystem<Input>();
    if (input->GetKeyPress(KEY_SPACE))
        animate_ = !animate_;

    // Move the camera, scale movement with time step
    MoveCamera(timeStep);

    // Animate objects' vertex data if enabled
    //if (animate_)
    //  AnimateObjects(timeStep);
}
