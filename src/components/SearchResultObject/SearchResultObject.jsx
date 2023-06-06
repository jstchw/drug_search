import React from "react";
import {Col, Container, Row, Badge} from "react-bootstrap";
import Histogram from "../Histogram/Histogram";
import './SearchResultObject.css'
import DrugDescription from "../DrugDescription/DrugDescription";


const SearchResultObject = ({ searchResults, searchTerm }) => {
    const { termCountDict, totalCount } = searchResults.result
    const [retrievedTermArr, setRetrievedTermArr] = React.useState(null)

    const processDrugGroups = (groups) => {
        return groups.map((group, index) => {
            let variant
            switch (group) {
                case 'approved':
                    variant = 'success'
                    break
                case 'investigational':
                    variant = 'warning'
                    break
                case 'illicit':
                    variant = 'danger'
                    break
                case 'vet_approved':
                    variant = 'info'
                    break
                default:
                    variant = 'secondary'
            }
            group = group.substring(0, 3).toUpperCase()
            return <Badge className={'mx-1'} key={index} bg={variant}>{group}</Badge>
        })

    }

    return (
        <div>
            <Container className="my-4">
                {retrievedTermArr && <h1>{retrievedTermArr[0]}</h1>}
                {retrievedTermArr && <h3>{processDrugGroups(retrievedTermArr[1])}</h3>}
            </Container>
            <Container className="">
                <Row className={'alert alert-secondary'}>
                    <span>The provided information is indicative and shouldn't be used for inference.</span>
                </Row>
                <Row>
                    <Col className="testing">
                        <DrugDescription drugName={searchTerm} onRetrieved={setRetrievedTermArr}/>
                    </Col>
                    <Col className={'testing'}>
                        {searchResults && <Histogram termCountDict={termCountDict} totalCount={totalCount} type={'all_groups'} />}
                    </Col>
                </Row>
            </Container>
        </div>
    );
}

export default SearchResultObject;