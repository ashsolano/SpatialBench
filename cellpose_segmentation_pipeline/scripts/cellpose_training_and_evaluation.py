from cellpose.io import logger_setup
from cellpose import io, metrics, models


# Initialse input training and test data 

train_images = []
train_labels = []
test_images = []
test_labels = []


# Initialise model and train on combined (raw and augmented) data
model_type = 'cyto2'

# Initialise Cellpose model
model = models.CellposeModel(gpu=True, model_type=model_type, net_avg=True, nchan=2)

logger_setup()

resample=True
net_average = True
channel_axis  = 0 

do_3D = False

#We are using our own normalisation
normalize = True

#Using the grey and channels
channels = [0,3] 

model = models.CellposeModel(gpu=True, model_type=model_type,net_avg=net_average,nchan=2)


model_save_path = "/path/to/save"
lr = 0.1
weight_decay=1e-5
n_epochs = 100
model_name = "model_name"
name = model_name

model.train(train_data=train_images,train_labels=train_labels,test_data=test_images,test_labels=test_labels,channels=channels,min_train_masks=0,
            normalize=normalize,save_path=model_save_path,learning_rate=lr,weight_decay=weight_decay,n_epochs=n_epochs,model_name=name)



# Evaluating the Cellpose model

model_precision ={}
model_ap={}

model_path = "path/to/model"

logger_setup()

model_type = model_path

model = models.CellposeModel(gpu=True, model_type=model_type,net_avg=True)

test_dir = "path/to/test_data"

output = io.load_train_test_data(test_dir,mask_filter="_seg.npy")
test_data, test_labels = output[:2]

channels = [0,3]

masks = model.eval(test_data, 
                    channels=channels,
                    cellprob_threshold=-5.5,
                    flow_threshold=0.95)[0]

ap = metrics.average_precision(test_labels, masks)[0]
print(ap[:,0].mean())

model_precision[model_name] = ap[:,0].mean()
model_ap[model_name] = ap