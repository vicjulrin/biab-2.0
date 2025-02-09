/**
 * BON in a Box - Script service
 * No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)
 *
 * The version of the OpenAPI document: 1.0.0
 * Contact: jean-michel.lord@mcgill.ca
 *
 * NOTE: This class is auto generated by OpenAPI Generator (https://openapi-generator.tech).
 * https://openapi-generator.tech
 * Do not edit the class manually.
 *
 */


import ApiClient from "../ApiClient";
import Info from '../model/Info';
import ScriptRunResult from '../model/ScriptRunResult';

/**
* Default service.
* @module api/DefaultApi
* @version 1.0.0
*/
export default class DefaultApi {

    /**
    * Constructs a new DefaultApi. 
    * @alias module:api/DefaultApi
    * @class
    * @param {module:ApiClient} [apiClient] Optional API client implementation to use,
    * default to {@link module:ApiClient#instance} if unspecified.
    */
    constructor(apiClient) {
        this.apiClient = apiClient || ApiClient.instance;
    }


    /**
     * Callback function to receive the result of the getPipelineInfo operation.
     * @callback module:api/DefaultApi~getPipelineInfoCallback
     * @param {String} error Error message, if any.
     * @param {module:model/Info} data The data returned by the service call.
     * @param {String} response The complete HTTP response.
     */

    /**
     * Get metadata about this pipeline
     * @param {String} descriptionPath Where to find the pipeline in ./pipeline folder.
     * @param {module:api/DefaultApi~getPipelineInfoCallback} callback The callback function, accepting three arguments: error, data, response
     * data is of type: {@link module:model/Info}
     */
    getPipelineInfo(descriptionPath, callback) {
      let postBody = null;
      // verify the required parameter 'descriptionPath' is set
      if (descriptionPath === undefined || descriptionPath === null) {
        throw new Error("Missing the required parameter 'descriptionPath' when calling getPipelineInfo");
      }

      let pathParams = {
        'descriptionPath': descriptionPath
      };
      let queryParams = {
      };
      let headerParams = {
      };
      let formParams = {
      };

      let authNames = [];
      let contentTypes = [];
      let accepts = ['application/json'];
      let returnType = Info;
      return this.apiClient.callApi(
        '/pipeline/{descriptionPath}/info', 'GET',
        pathParams, queryParams, headerParams, formParams, postBody,
        authNames, contentTypes, accepts, returnType, null, callback
      );
    }

    /**
     * Callback function to receive the result of the getPipelineOutputs operation.
     * @callback module:api/DefaultApi~getPipelineOutputsCallback
     * @param {String} error Error message, if any.
     * @param {Object.<String, {String: String}>} data The data returned by the service call.
     * @param {String} response The complete HTTP response.
     */

    /**
     * Get the output folders of the scripts composing this pipeline
     * @param {String} id Where to find the pipeline in ./script folder.
     * @param {module:api/DefaultApi~getPipelineOutputsCallback} callback The callback function, accepting three arguments: error, data, response
     * data is of type: {@link Object.<String, {String: String}>}
     */
    getPipelineOutputs(id, callback) {
      let postBody = null;
      // verify the required parameter 'id' is set
      if (id === undefined || id === null) {
        throw new Error("Missing the required parameter 'id' when calling getPipelineOutputs");
      }

      let pathParams = {
        'id': id
      };
      let queryParams = {
      };
      let headerParams = {
      };
      let formParams = {
      };

      let authNames = [];
      let contentTypes = [];
      let accepts = ['application/json'];
      let returnType = {'String': 'String'};
      return this.apiClient.callApi(
        '/pipeline/{id}/outputs', 'GET',
        pathParams, queryParams, headerParams, formParams, postBody,
        authNames, contentTypes, accepts, returnType, null, callback
      );
    }

    /**
     * Callback function to receive the result of the getScriptInfo operation.
     * @callback module:api/DefaultApi~getScriptInfoCallback
     * @param {String} error Error message, if any.
     * @param {module:model/Info} data The data returned by the service call.
     * @param {String} response The complete HTTP response.
     */

    /**
     * Get metadata about this script
     * @param {String} scriptPath Where to find the script in ./script folder.
     * @param {module:api/DefaultApi~getScriptInfoCallback} callback The callback function, accepting three arguments: error, data, response
     * data is of type: {@link module:model/Info}
     */
    getScriptInfo(scriptPath, callback) {
      let postBody = null;
      // verify the required parameter 'scriptPath' is set
      if (scriptPath === undefined || scriptPath === null) {
        throw new Error("Missing the required parameter 'scriptPath' when calling getScriptInfo");
      }

      let pathParams = {
        'scriptPath': scriptPath
      };
      let queryParams = {
      };
      let headerParams = {
      };
      let formParams = {
      };

      let authNames = [];
      let contentTypes = [];
      let accepts = ['application/json'];
      let returnType = Info;
      return this.apiClient.callApi(
        '/script/{scriptPath}/info', 'GET',
        pathParams, queryParams, headerParams, formParams, postBody,
        authNames, contentTypes, accepts, returnType, null, callback
      );
    }

    /**
     * Callback function to receive the result of the pipelineListGet operation.
     * @callback module:api/DefaultApi~pipelineListGetCallback
     * @param {String} error Error message, if any.
     * @param {Array.<String>} data The data returned by the service call.
     * @param {String} response The complete HTTP response.
     */

    /**
     * Get a list of available pipelines
     * @param {module:api/DefaultApi~pipelineListGetCallback} callback The callback function, accepting three arguments: error, data, response
     * data is of type: {@link Array.<String>}
     */
    pipelineListGet(callback) {
      let postBody = null;

      let pathParams = {
      };
      let queryParams = {
      };
      let headerParams = {
      };
      let formParams = {
      };

      let authNames = [];
      let contentTypes = [];
      let accepts = ['application/json'];
      let returnType = ['String'];
      return this.apiClient.callApi(
        '/pipeline/list', 'GET',
        pathParams, queryParams, headerParams, formParams, postBody,
        authNames, contentTypes, accepts, returnType, null, callback
      );
    }

    /**
     * Callback function to receive the result of the runPipeline operation.
     * @callback module:api/DefaultApi~runPipelineCallback
     * @param {String} error Error message, if any.
     * @param {String} data The data returned by the service call.
     * @param {String} response The complete HTTP response.
     */

    /**
     * Run this pipeline
     * @param {String} descriptionPath Where to find the script in ./script folder.
     * @param {Object} opts Optional parameters
     * @param {String} opts.body Content of input.json for this run
     * @param {module:api/DefaultApi~runPipelineCallback} callback The callback function, accepting three arguments: error, data, response
     * data is of type: {@link String}
     */
    runPipeline(descriptionPath, opts, callback) {
      opts = opts || {};
      let postBody = opts['body'];
      // verify the required parameter 'descriptionPath' is set
      if (descriptionPath === undefined || descriptionPath === null) {
        throw new Error("Missing the required parameter 'descriptionPath' when calling runPipeline");
      }

      let pathParams = {
        'descriptionPath': descriptionPath
      };
      let queryParams = {
      };
      let headerParams = {
      };
      let formParams = {
      };

      let authNames = [];
      let contentTypes = ['text/plain'];
      let accepts = ['text/plain'];
      let returnType = 'String';
      return this.apiClient.callApi(
        '/pipeline/{descriptionPath}/run', 'POST',
        pathParams, queryParams, headerParams, formParams, postBody,
        authNames, contentTypes, accepts, returnType, null, callback
      );
    }

    /**
     * Callback function to receive the result of the runScript operation.
     * @callback module:api/DefaultApi~runScriptCallback
     * @param {String} error Error message, if any.
     * @param {module:model/ScriptRunResult} data The data returned by the service call.
     * @param {String} response The complete HTTP response.
     */

    /**
     * Run this script
     * Run the script specified in the URL. Must include the extension.
     * @param {String} scriptPath Where to find the script in ./script folder
     * @param {Object} opts Optional parameters
     * @param {String} opts.body Content of input.json for this run
     * @param {module:api/DefaultApi~runScriptCallback} callback The callback function, accepting three arguments: error, data, response
     * data is of type: {@link module:model/ScriptRunResult}
     */
    runScript(scriptPath, opts, callback) {
      opts = opts || {};
      let postBody = opts['body'];
      // verify the required parameter 'scriptPath' is set
      if (scriptPath === undefined || scriptPath === null) {
        throw new Error("Missing the required parameter 'scriptPath' when calling runScript");
      }

      let pathParams = {
        'scriptPath': scriptPath
      };
      let queryParams = {
      };
      let headerParams = {
      };
      let formParams = {
      };

      let authNames = [];
      let contentTypes = ['text/plain'];
      let accepts = ['application/json'];
      let returnType = ScriptRunResult;
      return this.apiClient.callApi(
        '/script/{scriptPath}/run', 'POST',
        pathParams, queryParams, headerParams, formParams, postBody,
        authNames, contentTypes, accepts, returnType, null, callback
      );
    }

    /**
     * Callback function to receive the result of the scriptListGet operation.
     * @callback module:api/DefaultApi~scriptListGetCallback
     * @param {String} error Error message, if any.
     * @param {Array.<String>} data The data returned by the service call.
     * @param {String} response The complete HTTP response.
     */

    /**
     * Get a list of available scripts
     * @param {module:api/DefaultApi~scriptListGetCallback} callback The callback function, accepting three arguments: error, data, response
     * data is of type: {@link Array.<String>}
     */
    scriptListGet(callback) {
      let postBody = null;

      let pathParams = {
      };
      let queryParams = {
      };
      let headerParams = {
      };
      let formParams = {
      };

      let authNames = [];
      let contentTypes = [];
      let accepts = ['application/json'];
      let returnType = ['String'];
      return this.apiClient.callApi(
        '/script/list', 'GET',
        pathParams, queryParams, headerParams, formParams, postBody,
        authNames, contentTypes, accepts, returnType, null, callback
      );
    }

    /**
     * Callback function to receive the result of the stopPipeline operation.
     * @callback module:api/DefaultApi~stopPipelineCallback
     * @param {String} error Error message, if any.
     * @param data This operation does not return a value.
     * @param {String} response The complete HTTP response.
     */

    /**
     * Stop the specified pipeline run
     * @param {String} id Where to find the pipeline in ./script folder.
     * @param {module:api/DefaultApi~stopPipelineCallback} callback The callback function, accepting three arguments: error, data, response
     */
    stopPipeline(id, callback) {
      let postBody = null;
      // verify the required parameter 'id' is set
      if (id === undefined || id === null) {
        throw new Error("Missing the required parameter 'id' when calling stopPipeline");
      }

      let pathParams = {
        'id': id
      };
      let queryParams = {
      };
      let headerParams = {
      };
      let formParams = {
      };

      let authNames = [];
      let contentTypes = [];
      let accepts = [];
      let returnType = null;
      return this.apiClient.callApi(
        '/pipeline/{id}/stop', 'GET',
        pathParams, queryParams, headerParams, formParams, postBody,
        authNames, contentTypes, accepts, returnType, null, callback
      );
    }


}
